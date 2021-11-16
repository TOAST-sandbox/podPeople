#!/usr/bin/env bash

NSLOTS=${NSLOTS:-1}

set -eu

ROOT=$(dirname $(realpath $0))
export PATH=$ROOT:$PATH
export PATH=$ROOT/../sra-human-scrubber-1.1.2021-05-05/bin:$PATH
export PATH=$ROOT/../sra-human-scrubber-1.1.2021-05-05/scripts:$PATH

function check {
  set +e
  echo "Checking dependencies"
  echo "PATH: $PATH"
  for exe in bedtools bbsplitpairs.sh bowtie2 gcc samtools fasten_shuffle replaceReadsWithReference.pl scrub.sh; do
    which $exe 2>/dev/null
    if [ $? -gt 0 ]; then
      echo "Could not find $exe in PATH"
      exit
    fi
    $exe --help >/dev/null 2>&1
    if [ $? -gt 0 ]; then
      echo "ERROR with running $exe -h"
      exit 1
    fi
  done
  which aligns_to || exit 1
  exit 0
}
export -f check

ARGC=("$#")

# Check for dependencies if --check
if [[ "$1" =~ -+check ]]; then
  check
  exit 0
fi

if [[ "$ARGC" -lt 3 || "$1" == "--help" || "$1" == "-help" || "$1" == "-h" ]] ; then
   echo "
   ===========

   This simple script to replaces the sequence of human read contaminants from individuals (patients) with the corresponding sequence in a published human reference genome.
   
   ===========
   
   Usage: ./ReadCleanBot.general <in-dir> <out-dir> <path-to-ref>

   in-dir  Path to the directory containing the raw reads pairs in fastq format with filenames like 'name_S1_L001_R[12]_001.fastq'.
   
   out-dir  Path to the new output directory.
       
   path-to-ref  Path to the human reference genome assembly fasta (must be indexed for bowtie2)
   
   ===========
   
   " >&2 ;
   exit 0 ;
fi ;

### make output directories if needed
if [[ ! -d "$2" ]]; then
  mkdir $2;
  mkdir $2/02.fastq-interleaved
  mkdir $2/03.fastq-scrubbed
  mkdir $2/04.bam
  mkdir $2/05.fastq-replaceref
  mkdir $2/06.fastq-forSRA
  mkdir $2/07.summary
fi

#cd "${o}"

### wrap pipelin into a function 
function anonymize {
  r1=$1;
  outdir=$2
  ref=$3;    
  readid=$(basename $r1 | sed -e 's/_L[0-9]*.*//g' -e 's/_1.f.*//');

# step-1: Interleave R1 and R2

  r2=$(echo $r1 | sed -e 's/_R1_/_R2_/g' -e 's/_1.f/_2.f/');
  intlv="$outdir/02.fastq-interleaved/$readid.interleaved.fastq"
  #intlv=$(basename $r1 | sed 's/R1_001/interleaved/g');

  echo "R1          $r1"
  echo "R2          $r2"
  echo "interleaved $intlv"

  if [ ! -f $intlv ];then
    zcat $r1 $r2 | fasten_shuffle > $intlv;
  fi

# step-2: Human read scrubber and repair broken pairs

  scrb=$outdir/03.fastq-scrubbed/$(basename $intlv .interleaved.fastq).scrubbed.fastq;
  echo "scrubbed    $scrb"
  tmp=${scrb/.scrubbed.fastq/.tmp.fastq}
  single=${scrb/.scrubbed.fastq/.single.fastq}
  rtmp=${scrb/.scrubbed.fastq/.rtmp.fastq}
  rsingle=${scrb/.scrubbed.fastq/.rsingle.fastq}
  rmvd=${scrb/.scrubbed.fastq/.removed.fastq}

  if [ ! -f $scrb ];then
    scrub.sh -r $intlv;
    cat $intlv.clean | sed 's/\s/\//g' | sed 's/\:N\:0.*//g' > $tmp;
    cat $intlv.removed_spots | sed 's/\s/\//g' | sed 's/\:N\:0.*//g' > $rtmp;

    set -x
    bbsplitpairs.sh in=$tmp out=$scrb outs=$single fint
    bbsplitpairs.sh in=$rtmp out=$rmvd outs=$rsingle fint

  fi
  echo -e "\tScrubbed....[x]";

# step-3: Bowtie2 map to human reference

  bm=$outdir/04.bam/$(basename $intlv .fastq).t2t.bam;

  if [ ! -f $bm ];then
    bowtie2 \
    -X 1000 \
    -p "${NSLOTS}" -x $ref \
    --interleaved $rmvd | samtools view -b -F 12 > $bm.tmp;

    mv $bm.tmp $bm
  fi
  echo -e "\tMapped to reference....[x]";

# step-4: Replace patient sequence with corresponding T2T reference sequence

  rfq=$outdir/05.fastq-replaceref/$(basename $bm bam).fastq;

  if [ ! -f $rfq ];then
    replaceReadsWithReference.pl \
    $ref $bm 1> $rfq.tmp 2> $rfq.stderr.log;
    mv $rfq.tmp $rfq

    echo TODO
    exit 43
    
    # Not sure what this step is, because $rfq.clean doesn't exist
    #mv $outdir/05.fastq-replaceref/$rfq.clean $outdir/05.fastq-replaceref/$rfq.tmp

    bbsplitpairs.sh in=$outdir/05.fastq-replaceref/$rfq.tmp out=$outdir/05.fastq-replaceref/$rfq.clean outs=$outdir/05.fastq-replaceref/$rfq.single fint

  fi
  echo -e "\tRead sequences replaced....[x]";


# step-5: Combine, shuffle, separate into R1 and R2 for SRA submission

  f1=$(basename $r1);
  f2=$(basename $r2);

  if [ ! -f $outdir/06.fastq-forSRA/$f1 ];then
    cat $outdir/03.fastq-scrubbed/$scrb $outdir/05.fastq-replaceref/$rfq | sed 's/\sreplaced.*//g' | \
    fasten_randomize -n "${NSLOTS}" -p | \
    fasten_shuffle -n "${NSLOTS}" -d \
    -1 $outdir/06.fastq-forSRA/$f1 -2 $outdir/06.fastq-forSRA/$f2;
  fi
  echo -e "\tReads prepared for SRA....[x]";


# step-6: Summarize some count metrics

  sum=$(basename $intlv .fastq).stats.tsv;

  RAW=$(wc -l $outdir/02.fastq-interleaved/$intlv | awk '{print $1/8}');
  SCRpass=$(wc -l $outdir/03.fastq-scrubbed/$scrb | awk '{print $1/8}');
  T2Treplaced=$(wc -l $outdir/05.fastq-replaceref/$rfq | awk '{print $1/8}');
  TOTsra=$(wc -l $outdir/06.fastq-forSRA/$f1 | awk '{print $1/4}');
  FHS=$(echo "scale=4; $T2Treplaced/$TOTsra" | bc | awk '{printf "%.4f", $0}' );
  SCRfail=$(wc -l $outdir/03.fastq-scrubbed/$rmvd | awk '{print $1/8}');

  if [ "$RAW" -gt "$SCRpass" ];then
          T2F=$(echo "scale=4; $T2Treplaced/$SCRfail" | bc | awk '{printf "%.4f", $0}' );
  else
          T2F=0;
  fi

  echo -e "ReadID\tRaw_pairs\tScrub_pass\tScrub_remove\tReplaced_pairs\tReplaced_frac\tTotal_pairs_final\tHuman_frac_final" > ./07.summary/$sum;
  echo -e "$readid\t$RAW\t$SCRpass\t$SCRfail\t$T2Treplaced\t$T2F\t$TOTsra\t$FHS" >> ./07.summary/$sum;


echo -e "\tPipeline finished!\n";

}
export -f anonymize
####

echo -e "\nInput directory is:" $1;
echo -e "Output directory is:" $2;
echo -e "\nAnonymizing host reads..."; 

#for i in $(find $1/*R1*fastq); do
R1s=$(\ls -f1 $1/*_R1*.fastq.gz $1/*_1.fastq.gz || true)
for i in $R1s; do
  anonymize $i $2 $3;

done

cat $2/07.summary/*tsv | head -1 > $2/Summary.tsv;
cat $2/07.summary/*tsv | grep -v "^ReadID" | sort >> $2/Summary.tsv;

