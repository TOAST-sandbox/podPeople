#!/usr/bin/env bash

NSLOTS=${NSLOTS:-1}

set -eu

ROOT=$(dirname $(realpath $0))
export PATH=$ROOT:$PATH
export PATH=$ROOT/../sra-human-scrubber-1.1.2021-05-05/bin:$PATH
export PATH=$ROOT/../sra-human-scrubber-1.1.2021-05-05/scripts:$PATH

function check {
  # Handling errors ourselves within this function
  set +e

  echo "Checking dependencies"
  echo "PATH: $PATH"
  dependencies=(
    bedtools
    bbsplitpairs.sh
    bowtie2
    gcc
    samtools
    fasten_shuffle
    replaceReadsWithReference.pl
    scrub.sh
    )
  for exe in "${dependencies[@]}"; do
    which $exe 2> /dev/null
    if [ $? -gt 0 ]; then
      echo "ERROR: Could not find $exe in PATH" >&2
      exit 1
    fi
    $exe --help >/dev/null 2>&1
    if [ $? -gt 0 ]; then
      echo "ERROR: non-zero exit status running \`$exe --help\`" >&2
      exit 1
    fi
  done
  which aligns_to || exit 1
  exit 0
}
export -f check

ARGC=("$#")

if [[ "$ARGC" -lt 3 || "$1" == "--help" || "$1" == "-help" || "$1" == "-h" || "$1" == "" ]]; then
   echo "
===========

Replaces the sequence of human read contaminants from individuals (patients)
with the corresponding sequence in a published human reference genome.

===========

Usage: ./ReadCleanBot.general [--check] <in-dir> <out-dir> <ref-file>

Required:
  <in-dir>      Path to the directory containing the raw reads pairs in FastQ
                format with filenames like 'name_S1_L001_R[12]_001.fastq'
                or 'name_[12].fastq'.

  <out-dir>     Path to the new output directory.

  <ref-file>    Path to the human reference genome assembly FastA (must be
                indexed for bowtie2).

Optional:
   --check      Check dependencies before running the algorithm. Exit code is
                0 for no error; 1 for missing dependency.
   -h | --help  Show this help message and exit.

===========

   " >&2
   exit 0
fi

# Check for dependencies if --check
if [[ "$1" =~ -+check ]]; then
  check
  exit 0
fi


### make output directories if needed
outdirs=(
  "$2"
  "$2/02.fastq-interleaved"
  "$2/03.fastq-scrubbed"
  "$2/04.bam"
  "$2/05.fastq-replaceref"
  "$2/06.fastq-forSRA"
  "$2/07.summary"
  )
for outdir in "${outdirs[@]}"; do
  if [[ ! -d "$outdir" ]]; then
    mkdir $outdir
  fi
done

### wrap pipeline into a function
function anonymize {
  r1=$1
  outdir=$2
  ref=$3
  readid=$(basename $r1 | sed -e 's/_L[0-9]*.*//g' -e 's/_1.f.*//')

# step-1: Interleave R1 and R2

  r2=$(echo $r1 | sed -e 's/_R1_/_R2_/g' -e 's/_1.f/_2.f/')
  if [ ! -f "${r2}" ]; then
    echo "ERROR: ${r2} absent. Unable to pair with ${r1}" >&2
    exit 1
  fi

  intlv="$outdir/02.fastq-interleaved/$readid.interleaved.fastq"

  echo "R1          $r1"
  echo "R2          $r2"
  echo "interleaved $intlv"

  if [ ! -f $intlv ]; then
    zcat $r1 $r2 | fasten_shuffle > $intlv
  fi

# step-2: Human read scrubber and repair broken pairs

  scrb=$outdir/03.fastq-scrubbed/$(basename $intlv .interleaved.fastq).scrubbed.fastq
  tmp=${scrb/.scrubbed.fastq/.tmp.fastq}
  single=${scrb/.scrubbed.fastq/.single.fastq}
  rtmp=${scrb/.scrubbed.fastq/.rtmp.fastq}
  rsingle=${scrb/.scrubbed.fastq/.rsingle.fastq}
  rmvd=${scrb/.scrubbed.fastq/.removed.fastq}

  if [ ! -f $scrb ]; then
    scrub.sh -r $intlv
    cat $intlv.clean | sed 's/\s/\//g' | sed 's/\:N\:0.*//g' > $tmp
    cat $intlv.removed_spots | sed 's/\s/\//g' | sed 's/\:N\:0.*//g' > $rtmp

    bbsplitpairs.sh in=$tmp out=$scrb outs=$single fint
    bbsplitpairs.sh in=$rtmp out=$rmvd outs=$rsingle fint
  fi
  echo -e "\tScrubbed....[x]"

# step-3: Bowtie2 map to human reference

  bm=$outdir/04.bam/$(basename $intlv .fastq).t2t.bam

  if [ ! -f $bm ]; then
    set -x
    bowtie2 \
     -X 1000 \
     -p "${NSLOTS}" \
     -x $ref \
     --interleaved $scrb | \
     samtools view -b | \
     samtools sort -n - > $bm.tmp

    mv $bm.tmp $bm
    set +x
  fi
  echo -e "\tMapped to reference....[x]"

# step-4: Replace patient sequence with corresponding T2T reference sequence

  rfq=$outdir/05.fastq-replaceref/$(basename $bm .bam).fastq

  if [ ! -f $rfq ]; then
    set -x
    replaceReadsWithReference.pl \
     $ref $bm 1> $rfq.tmp 2> $rfq.stderr.log

    mv $rfq.tmp $rfq
    #bbsplitpairs.sh in=$rfq.tmp out=$rfq.clean outs=$rfq.single fint
    set +x
  fi
  echo -e "\tRead sequences replaced....[x]"


# step-5: Combine, shuffle, separate into R1 and R2 for SRA submission

  f1=$outdir/06.fastq-forSRA/$(basename $r1)
  f2=$outdir/06.fastq-forSRA/$(basename $r2)

  if [ ! -f $f1 ]; then
    cat $scrb $rfq | \
     sed 's/\sreplaced.*//g' | \
     fasten_randomize -n "${NSLOTS}" -p | \
     fasten_shuffle -n "${NSLOTS}" -d -1 $f1.tmp -2 $f2.tmp

    mv $f1.tmp $f1
    mv $f2.tmp $f2
  fi
  echo -e "\tReads prepared for SRA....[x]"

# step-6: Summarize some count metrics

  sum=$outdir/07.summary/$(basename $intlv .fastq).stats.tsv

  RAW=$(wc -l $intlv | awk '{print $1/8}')
  SCRpass=$(wc -l $scrb | awk '{print $1/8}')
  T2Treplaced=$(wc -l $rfq | awk '{print $1/8}')
  TOTsra=$(wc -l $f1 | awk '{print $1/4}')
  FHS=$(echo "scale=4; $T2Treplaced/$TOTsra" | bc | awk '{printf "%.4f", $0}')
  SCRfail=$(wc -l $rmvd | awk '{print $1/8}')

  if [ "$RAW" -gt "$SCRpass" ]; then
    T2F=$(echo "scale=4; $T2Treplaced/$SCRfail" | bc | awk '{printf "%.4f", $0}')
  else
    T2F=0
  fi

  echo -e "ReadID\tRaw_pairs\tScrub_pass\tScrub_remove\tReplaced_pairs\tReplaced_frac\tTotal_pairs_final\tHuman_frac_final" > $sum
  echo -e "$readid\t$RAW\t$SCRpass\t$SCRfail\t$T2Treplaced\t$T2F\t$TOTsra\t$FHS" >> $sum

echo -e "\tPipeline finished!\n"

}
export -f anonymize
####

echo -e "\nInput directory is:" $1
echo -e "Output directory is:" $2
echo -e "\nAnonymizing host reads..."

R1s=$(\ls -f1 $1/*_R1*.fastq.gz $1/*_1.fastq.gz || true)
for i in $R1s; do
  anonymize $i $2 $3
done

cat $2/07.summary/*.tsv | head -1 > $2/Summary.tsv
cat $2/07.summary/*.tsv | grep -v "^ReadID" | sort >> $2/Summary.tsv
