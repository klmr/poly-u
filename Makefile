SHELL := $(shell which bash)
ref_dir := data/reference
index_dir := data/index
reference := ${ref_dir}/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa
viral-reference := ${ref_dir}/orsay-virus.fa
infected-reference := ${ref_dir}/$(basename $(notdir ${reference}))-orv.fa
annotation := ${ref_dir}/Caenorhabditis_elegans.WBcel235.84.gtf
gene-annotation := ${ref_dir}/Caenorhabditis_elegans.WBcel235.84.genes.gtf
viral-annotation := ${ref_dir}/orsay-virus.gtf
infected-gene-annotation := ${ref_dir}/Caenorhabditis_elegans.WBcel235.84.genes-orv.gtf
index := ${index_dir}/Caenorhabditis_elegans/Genome
viral-index := ${index_dir}/orsay-virus/Genome
infected-index := ${index_dir}/Caenorhabditis_elegans-orv/Genome

fastq = $(foreach r,R5 R3,raw/$(shell grep --only-matching c_elegans_.. <<< "$1")/fastq/$(notdir $(subst _R5,_$r,$1)))

bsub := scripts/bsub -K -q research-rh7

.PHONY: reference
reference: ${reference}

${reference}:
	mkdir -p "$(dir $@)"
	wget 'ftp://ftp.ensembl.org/pub/release-84/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz' -O $@.gz
	gunzip $@.gz

.PHONY: viral-reference
viral-reference: ${viral-reference}

${viral-reference}: raw/reference/$(notdir ${viral-reference})
	mkdir -p "$(dir $@)"
	cp "$<" "$@"

.PHONY: infected-reference
infected-reference: ${infected-reference}

${infected-reference}: ${reference} ${viral-reference}
	cat $+ > '$@'

.PHONY: annotation
annotation: ${annotation}

${annotation}:
	mkdir -p "$(dir $@)"
	wget 'ftp://ftp.ensembl.org/pub/release-84/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.84.gtf.gz' -O $@.gz
	gunzip $@.gz

.PHONY: viral-annotation
viral-annotation: ${viral-annotation}

${viral-annotation}: raw/reference/orsay.genome.ape
	mkdir -p "$(dir $@)"
	./scripts/genbank-to-gtf "$<" "$@"

# Alignment

.PHONY: index
index: ${index}

${index}: ${reference}
	mkdir -p "$(dir $@)"
	${bsub} -n 12 "STAR --runThreadN 12 --runMode genomeGenerate \
		--genomeDir '$(dir $@)' --genomeFastaFiles '$<'"
	rm Log.out

.PHONY: viral-index
viral-index: ${viral-index}

${viral-index}: ${viral-reference}
	mkdir -p "$(dir $@)"
	STAR --runMode genomeGenerate --genomeSAindexNbases 3 \
		--genomeDir '$(dir $@)' --genomeFastaFiles '$<'
	rm Log.out

.PHONY: infected-index
infected-index: ${infected-index}

${infected-index}: ${infected-reference}
	mkdir -p "$(dir $@)"
	${bsub} -n 12 "STAR --runThreadN 12 --runMode genomeGenerate \
		--genomeDir '$(dir $@)' --genomeFastaFiles '$<'"
	rm Log.out

raw-reads = $(shell ls raw/c_elegans_??/fastq/*_R?.fastq.gz)
trimmed-reads = $(filter-out %_R3.fastq.gz,$(subst /fastq/,/,$(subst raw/,data/trimmed/,${raw-reads})))

.PHONY: trimmed-reads
trimmed-reads: ${trimmed-reads}

.SECONDEXPANSION:

data/trimmed/%_R5.fastq.gz: $$(call fastq,$$@)
	mkdir -p data/qc
	mkdir -p "$(dir $@)"
	${bsub} "cutadapt -a TGGAATTCTCGG -A TGGAATTCTCGG -G GTTCAGAGTTCTACAGTCCGACGATC \
		--minimum-length 5 -o '$@' -p '$(subst _R5,_R3,$@)' $^"
	fastqc -o data/qc '$@'
	fastqc -o data/qc '$(subst _R5,_R3,$@)'

mapped-reads = $(subst /trimmed/,/mapped/,$(filter-out %_R3.fastq.gz,${trimmed-reads:_R5.fastq.gz=.bam}))

.PHONY: mapped-reads
mapped-reads: ${mapped-reads}

data/mapped/%.bam: data/trimmed/%_R5.fastq.gz ${infected-index}
	mkdir -p "$(dir $@)"
	${bsub} -n 12 -M24000 -R'select[mem>24000]' -R'rusage[mem=24000]' \
		"STAR --runThreadN 12 --genomeDir '$(dir ${infected-index})' \
		--runMode alignReads --alignEndsType EndToEnd \
		--scoreInsOpen -10000 --scoreDelOpen -10000 \
		--outFilterMismatchNmax 0 --outFilterMultimapNmax 1 \
		--readFilesIn $< --readFilesCommand 'gunzip -c' \
		--outSAMtype BAM Unsorted --outFileNamePrefix '$(basename $@)'"
	mv "$(basename $@)Aligned.out.bam" "$(basename $@).bam"

.PHONY: qc-report
## Generate aggregate report from individual tool/QC outputs
qc-report: data/qc/multiqc_report.html

data/qc/multiqc_report.html: ${trimmed-reads} ${mapped-reads}
	multiqc --force --outdir data/qc \
		data/trimmed data/qc data/mapped

data/softclip-mapped/%.bam: data/trimmed/%_R5.fastq.gz ${infected-index}
	mkdir -p "$(dir $@)"
	${bsub} -n 12 -M24000 -R'select[mem>24000]' -R'rusage[mem=24000]' \
		"STAR --runThreadN 12 --genomeDir '$(dir ${infected-index})' \
		--runMode alignReads --alignEndsType Local \
		--outFilterMultimapNmax 1 \
		--readFilesIn $< --readFilesCommand 'gunzip -c' \
		--outSAMtype BAM Unsorted --outFileNamePrefix '$(basename $@)'"
	mv "$(basename $@)Aligned.out.bam" "$(basename $@).bam"

softclip-indexed = $(subst /mapped/,/softclip-mapped/,${mapped-reads:.bam=-sorted.bam.bai})
.SECONDARY: ${softclip-indexed}

data/softclip-mapped/%-sorted.bam.bai: data/softclip-mapped/%.bam
	samtools sort -o "$(basename $@)" "$<"
	samtools index "$(basename $@)"

# The mapped data with soft-clipping contains many spurious matches to the
# viral RNA. We filter this by taking all putative viral RNA hits and verify
# that their 3p tails actually align well to the 3p end of the viral RNA
# reference.

viral-correction = $(subst /mapped/,/viral/,${mapped-reads:.bam=.tsv})

.PHONY: viral-correction
viral-correction: ${viral-correction}

.SECONDARY: ${viral-correction}

data/viral/%.tsv: data/softclip-mapped/%-sorted.bam.bai data/trimmed/%_R3.fastq.gz ${viral-reference}
	mkdir -p "$(dir $@)"
	samtools view '${<:.bai=}' ORV-RNA1 ORV-RNA2 | cut -f1,3 > '$(basename $@).id'
	${bsub} "./scripts/verify-3p-mapping --fastq '$(word 2,$^)' \
		--reference '$(lastword $^)' '$(basename $@).id' > '$@'"
	rm -f '$(basename $@).id'

${gene-annotation}: ${annotation}
	awk '($$3 == "gene") {print $$0}' '$<' > '$@'

${infected-gene-annotation}: ${gene-annotation} ${viral-annotation}
	cat $+ > '$@'

find-genes := $(subst /mapped/,/genes/,${mapped-reads:.bam=.tsv})

.PHONY: find-genes
find-genes: ${find-genes}

data/genes/%.tsv: data/mapped/%.bam ${infected-gene-annotation}
	mkdir -p "$(dir $@)"
	${bsub} "./scripts/find-mapped-genes '$<' '$@' '${infected-gene-annotation}'"

aligned-taginfo := $(subst /genes/,/taginfo/aligned/,${find-genes})

.PHONY: aligned-taginfo
aligned-taginfo: ${aligned-taginfo}

# trimmed-fastq_r3 = $(subst /genes/,/trimmed/,${1:.tsv=_R3.fastq.gz})

# data/taginfo/aligned/%.tsv: data/genes/%.tsv ${infected-reference} ${infected-gene-annotation}
# 	mkdir -p "$(dir $@)"
# 	${bsub} -n 16 -M 24000 -R 'span[hosts=1] select[mem>24000] rusage[mem=24000]' \
# 		"./scripts/3p-align --reference '${infected-reference}' \
# 		--annotation '${infected-gene-annotation}' \
# 		--genes '$<' --ncores 16 '$(call trimmed-fastq_r3,$<)' > '$@'"

taginfo = $(subst /genes/,/taginfo/,${find-genes})

.PHONY: taginfo
taginfo: ${taginfo}

data/taginfo/%.tsv: data/genes/%.tsv data/trimmed/%_R3.fastq.gz data/viral/%.tsv
	mkdir -p "$(dir $@)"
	${bsub} -M 16000 -R'select[mem>16000] rusage[mem=16000]' \
		"./scripts/merge-taginfo \
		--genes '$<' --reads '$(word 2,$^)' --viral '$(lastword $^)' > '$@'"

.PHONY: merged-taginfo
merged-taginfo: data/taginfo/all-tailinfo.tsv

data/taginfo/all-taginfo.tsv: ${taginfo}
	head -n 1 "$<" > "$@"
	tail -q -n +2 $+ >> "$@"

data/taginfo/all-tailinfo.tsv: data/taginfo/all-taginfo.tsv
	./scripts/find-tail-modifications '$<' > '$@'

data/taginfo/poly-a-summary.tsv: data/taginfo/all-tailinfo.tsv
	./scripts/summarize-poly-a-length '$<' > '$@'

data/taginfo/tail-mod-summary.tsv: data/taginfo/all-tailinfo.tsv
	./scripts/summarize-tailinfo '$<' > '$@'

.PHONY: summaries
summaries: data/taginfo/poly-a-summary.tsv data/taginfo/tail-mod-summary.tsv

plots = data/plots/uninfected-poly-a-lengths-density.pdf \
		data/plots/infected-poly-a-lengths-density.pdf \
		data/plots/uninfected-poly-a-lengths-gene-sets-density.pdf \
		data/plots/infected-poly-a-lengths-gene-sets-density.pdf \
		data/plots/uninfected-a-u-scatter.pdf \
		data/plots/infected-a-u-scatter.pdf \
		data/plots/viral-poly-u-lengths.pdf

.PHONY: plots
plots: ${plots}

data/plots/%-poly-a-lengths-density.pdf: data/taginfo/all-taginfo.tsv
	mkdir -p "$(dir $@)"
	./scripts/plot-global-poly-a-lengths --plot '$@' --treatment '$*' '$<'

data/plots/%-poly-a-lengths-gene-sets-density.pdf: data/taginfo/all-taginfo.tsv
	mkdir -p "$(dir $@)"
	./scripts/plot-poly-a-gene-sets --plot '$@' --treatment '$*' \
		--germline 'raw/Germline enriched genes from Reinke et al 2004 supp fig1.tsv' \
		--infection-response 'raw/Genes upregulated in inf rde-1 vs inf N2 from Sarkies et al 2012 supp table1.tsv' \
		'$<'

data/plots/%-a-u-scatter.pdf: data/taginfo/all-tailinfo.tsv
	mkdir -p "$(dir $@)"
	./scripts/plot-poly-a-u-correlations --plot '$@' --treatment '$*' '$<'

data/plots/viral-poly-u-lengths.pdf: data/taginfo/tail-mod-summary.tsv
	mkdir -p "$(dir $@)"
	./scripts/plot-viral-poly-u-lengths --plot '$@' '$<'
