data_dir := data
ref_dir := ${data_dir}/reference
index_dir := ${data_dir}/index
reference := ${ref_dir}/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa
annotation := ${ref_dir}/Caenorhabditis_elegans.WBcel235.84.gtf
gene-annotation := ${ref_dir}/Caenorhabditis_elegans.WBcel235.84.genes.gtf
index := ${index_dir}/Caenorhabditis_elegans/Genome

define find-fastq=
$(foreach r,R5 R3,raw/$(shell grep --only-matching c_elegans_.. <<< "$1")/fastq/$(basename $(notdir $1))_$r.fastq.gz)
endef

bsub := scripts/bsub -K -q research-rh7

.PHONY: reference
reference: ${reference}

${reference}:
	mkdir -p "$(dir $@)"
	wget 'ftp://ftp.ensembl.org/pub/release-84/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz' -O $@.gz
	gunzip $@.gz

.PHONY: annotation
annotation: ${annotation}

${annotation}:
	mkdir -p "$(dir $@)"
	wget 'ftp://ftp.ensembl.org/pub/release-84/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.84.gtf.gz' -O $@.gz
	gunzip $@.gz

# Alignment

.PHONY: index
index: ${index}

${index}: ${reference}
	mkdir -p "$(dir $@)"
	${bsub} -n 12 "STAR --runThreadN 4 --runMode genomeGenerate \
		--genomeDir '$(dir $@)' --genomeFastaFiles '$<'"

mapped-reads := $(foreach f,$(shell ls raw/c_elegans_*/fastq/*_R5.fastq.gz),${data_dir}/mapped/$(subst raw/,,$(subst fastq/,,$(subst _R5.fastq.gz,,$f))).bam)

.PHONY: mapped-reads
mapped-reads: ${mapped-reads}

.SECONDEXPANSION:

${data_dir}/mapped-paired-end/%.bam: $$(call find-fastq,%) ${index} ${annotation}
	mkdir -p "$(dir $@)"
	${bsub} -n 6 -M24000 -R'select[mem>24000]' -R'rusage[mem=24000]' \
		"STAR --runThreadN 4 --genomeDir '$(dir ${index})' \
		--runMode alignReads --alignEndsType Local \
		--sjdbGTFfile '${annotation}' \
		--readFilesIn $(call find-fastq,$@) --readFilesCommand 'gunzip -c' \
		--outSAMtype BAM Unsorted --outFileNamePrefix '$(basename $@).'"
	mv "$(basename $@).Aligned.out.bam" "$(basename $@).bam"
 
find-fastq_r5=raw/$(shell grep --only-matching c_elegans_.. <<< "$1")/fastq/$(basename $(notdir $1))_R5.fastq.gz

${data_dir}/mapped/%.bam: $$(call find-fastq_r5,%) ${index} ${annotation}
	mkdir -p "$(dir $@)"
	${bsub} -n 6 -M24000 -R'select[mem>24000]' -R'rusage[mem=24000]' \
		"STAR --runThreadN 4 --genomeDir '$(dir ${index})' \
		--runMode alignReads --alignEndsType Local \
		--sjdbGTFfile '${annotation}' \
		--readFilesIn $(call find-fastq_r5,$@) --readFilesCommand 'gunzip -c' \
		--outSAMtype BAM Unsorted --outFileNamePrefix '$(basename $@).'"
	mv "$(basename $@).Aligned.out.bam" "$(basename $@).bam"

${gene-annotation}: ${annotation}
	awk '($$3 == "gene") {print $$0}' '$<' > '$@'

find-genes := $(subst /mapped/,/genes/,${mapped-reads:.bam=.tsv})

.PHONY: find-genes
find-genes: ${find-genes}

${data_dir}/genes/%.tsv: ${data_dir}/mapped/%.bam ${gene-annotation}
	mkdir -p "$(dir $@)"
	${bsub} "./scripts/find-mapped-genes '$<' '$@' '${gene-annotation}'"
