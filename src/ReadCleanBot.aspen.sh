#!/bin/bash



# Setup new workspace and bring over initial raw read files
outpath=/scicomp/groups-pure/Projects/TOAST/bad.genomes/renamed.bad.genomes/cGul-another-round
cd "${outpath}"
# mkdir 01.fastq-renamed
# cp -r /scicomp/groups-pure/Projects/TOAST/bad.genomes/renamed.bad.genomes/mrwReformatted/01.fastq-renamed/*.fastq \
#  01.fastq-renamed/

# # Submit each read pair to run independently on the cluster
# mkdir 02.fastq-interleaved/ 03.fastq-scrubbed/ 04.bam-t2t/ 05.fastq-replaceref/ 06.fastq-forSRA/ 07.summary/
# mkdir log
for r1 in 01.fastq-renamed/bad*R1*.fastq; do
  b=$(basename "${r1}" | cut -d _ -f 1)
  qsub \
   -q all.q -q short.q -pe smp 6-12 \
   -N "${b}" \
   -e "${outpath}"/log \
   -o "${outpath}"/log \
   -v r1="${r1}",o="${outpath}" \
   "${outpath}"/ReadCleanBot.individual-sample.sh
done
