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
		--runMode alignReads --alignEndsType Local \
		--outFilterMultimapNmax 1 \
		--readFilesIn $< --readFilesCommand 'gunzip -c' \
		--outSAMtype BAM Unsorted --outFileNamePrefix '$(basename $@)'"
	mv "$(basename $@)Aligned.out.bam" "$(basename $@).bam"

.PHONY: qc-report
## Generate aggregate report from individual tool/QC outputs
qc-report: data/qc/multiqc_report.html

data/qc/multiqc_report.html: ${trimmed-reads} ${mapped-reads}
	multiqc --force --outdir data/qc \
		data/trimmed data/qc data/mapped

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

trimmed-fastq_r3 = $(subst /genes/,/trimmed/,${1:.tsv=_R3.fastq.gz})

data/taginfo/aligned/%.tsv: data/genes/%.tsv ${infected-reference} ${infected-gene-annotation}
	mkdir -p "$(dir $@)"
	${bsub} -n 16 -M 24000 -R 'span[hosts=1] select[mem>24000] rusage[mem=24000]' \
		"./scripts/3p-align --reference '${infected-reference}' \
		--annotation '${infected-gene-annotation}' \
		--genes '$<' --ncores 16 '$(call trimmed-fastq_r3,$<)' > '$@'"

taginfo = $(subst /genes/,/taginfo/,${find-genes})

.PHONY: taginfo
taginfo: ${taginfo}

find-reads_3p = raw/$(shell grep --only-matching c_elegans_.. <<< "$1")/fastq/$(notdir $1)_R3.fastq.gz

data/taginfo/%.tsv: data/genes/%.tsv $$(call find-reads_3p,%)
	mkdir -p "$(dir $@)"
	${bsub} "./scripts/merge-taginfo \
		--genes '$(firstword $+)' --reads '$(lastword $+)' > '$@'"

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
