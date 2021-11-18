.DELETE_ON_ERROR:

.LOW_RESOLUTION_TIME:

.DEFAULT: install db/t2t-chm13.fasta

install: ngs-tools/tools/tax/Makefile
	cd ngs-tools/tools/tax && ./quickbuild.sh
	@echo "remember: export PATH=$$(pwd -P)/ngs-tools/tools/tax/bin"

ngs-tools/tools/tax/Makefile:
	git clone https://github.com/ncbi/ngs-tools.git --branch tax

sra-human-scrubber-1.1.2021-05-05/data/human_filter.db:
	wget https://github.com/ncbi/sra-human-scrubber/archive/refs/tags/1.1.2021-05-05.tar.gz 
	tar zxvf 1.1.2021-05-05.tar.gz
	cd sra-human-scrubber-1.1.2021-05-05 && bash init_db.sh

# List of chromosomes obtained by:
#   esearch -db assembly -query GCA_009914755.3 | elink -target nuccore | esummary | xtract -pattern DocumentSummary -element Caption | perl -lane 'print "db/$F[0].chrom.fasta"' | tr '\n' ' '; echo;
db/t2t-chm13.fasta: db/CP068255.chrom.fasta db/CP068256.chrom.fasta db/CP068257.chrom.fasta db/CP068258.chrom.fasta db/CP068259.chrom.fasta db/CP068260.chrom.fasta db/CP068261.chrom.fasta db/CP068262.chrom.fasta db/CP068263.chrom.fasta db/CP068264.chrom.fasta db/CP068265.chrom.fasta db/CP068266.chrom.fasta db/CP068267.chrom.fasta db/CP068268.chrom.fasta db/CP068269.chrom.fasta db/CP068270.chrom.fasta db/CP068271.chrom.fasta db/CP068272.chrom.fasta db/CP068273.chrom.fasta db/CP068274.chrom.fasta db/CP068275.chrom.fasta db/CP068276.chrom.fasta db/CP068277.chrom.fasta db/CP068254.chrom.fasta
	cat $^ > $@
	grep '>' $^
	# We expect 22 autosomal + 1 X + 1 mitochondrial chromosome
	[[ "$$(grep -c '>' $@)" == 24 ]]
	# truncate the separate files but keep the filenames
	for i in $^; do echo -n > $$i; done;

db/t2t-chm13.fasta.1.bt2: db/t2t-chm13.fasta
	bowtie2-build $< $<

db/t2t-chm13.fasta.fai: db/t2t-chm13.fasta
	samtools faidx $<

db/%.chrom.fasta:
	mkdir -p $$(dirname $@)
	esearch -db nuccore -query $$(basename $@ .chrom.fasta) | efetch -format fasta > $@
	# Must have at least 2 lines in the fasta file
	# Use `head` so that it doesn't load the whole huge file
	if [[ "$$(head $@ | wc -l)" -lt 2 ]]; then echo "fasta file was fewer than 2 lines"; exit 1; fi;

