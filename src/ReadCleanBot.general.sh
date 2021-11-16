#!/usr/bin/env bash


if [[ "$1" == "" ||"$2" == "" ||"$3" == "" || "$1" == "-help" || "$1" == "-h" ]] ; then
   echo "
   ===========

   This simple script to replaces the sequence of human read contaminants from individuals (patients) with the corresponding sequence in a published human reference genome.
   
   ===========
   
   Usage: ./ReadCleanBot.general <in-dir> <out-dir> <path-to-ref>

   in-dir	Path to the directory containing the raw reads pairs in fastq format with filenames like 'name_S1_L001_R[12]_001.fastq'.
   
   out-dir	Path to the new output directory.
   		
   path-to-ref	Path to the human reference genome assembly fasta (must be indexed for bowtie2)
   
   ===========
   
   Dependenecies in your $PATH:	BBMap, BEDTools, bowtie2, samtools, fasten, replaceReadsWithReference.pl (lskScripts), scrub.sh (NCBI)  
   								(on biolinux - 'ml  BBMap/38.90 BEDTools/2.27.1 bowtie2/2.3.5.1 gcc/4.9.3 samtools/1.10')
     
   ===========
   " >&2 ;
   exit 1 ;
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
	ref=$3;		
	readid=$(basename $r1 | sed 's/_L001.*//g');

# step-1: Interleave R1 and R2

	r2=$(echo $r1 | sed 's/_R1_/_R2_/g');
	intlv=$(basename $r1 | sed 's/R1_001/interleaved/g');

	echo -e "$r1\t$r2\t$intlv";

	if [ ! -f $2/02.fastq-interleaved/$intlv ];then
		fasten_shuffle -n "${NSLOTS}" -1 $r1 -2 $r2 > \
		$2/02.fastq-interleaved/$intlv;
	fi

# step-2: Human read scrubber and repair broken pairs

	scrb=$(basename $intlv .fastq)_scrubbed.fastq;
	tmp=$(echo $scrb | sed 's/scrubbed/tmp/g');
	single=$(echo $scrb |sed 's/scrubbed/single/g');
    rtmp=$(echo $scrb | sed 's/scrubbed/rtmp/g');
    rsingle=$(echo $scrb |sed 's/scrubbed/rsingle/g');
    rmvd=$(echo $scrb |sed 's/scrubbed/removed/g');

	if [ ! -f $2/03.fastq-scrubbed/$scrb ];then
		scrub.sh -r $2/02.fastq-interleaved/$intlv;
		cat $2/02.fastq-interleaved/$intlv.clean | sed 's/\s/\//g' | sed 's/\:N\:0.*//g' > $2/03.fastq-scrubbed/$tmp;
		cat $2/02.fastq-interleaved/$intlv.removed_spots | sed 's/\s/\//g' | sed 's/\:N\:0.*//g' > $2/03.fastq-scrubbed/$rtmp;

		bbsplitpairs.sh in=$2/03.fastq-scrubbed/$tmp out=$2/03.fastq-scrubbed/$scrb outs=$2/03.fastq-scrubbed/$single fint
        bbsplitpairs.sh in=$2/03.fastq-scrubbed/$rtmp out=$2/03.fastq-scrubbed/$rmvd outs=$2/03.fastq-scrubbed/$rsingle fint

	fi
	echo -e "\tScrubbed....[x]";

# step-3: Bowtie2 map to human reference

	bm=$(basename $intlv .fastq).t2t.bam;

	if [ ! -f $2/04.bam/$bm ];then
		bowtie2 \
		-X 1000 \
		-p "${NSLOTS}" -x $ref \
                --interleaved $2/03.fastq-scrubbed/$rmvd | samtools view -b -F 12 > $2/04.bam/$bm;

	fi
	echo -e "\tMapped to reference....[x]";

# step-4: Replace patient sequence with corresponding T2T reference sequence

	rfq=$(basename $bm bam)fastq;

	if [ ! -f $2/05.fastq-replaceref/$rfq ];then
		replaceReadsWithReference.pl \
		$ref $2/04.bam/$bm 1> $2/05.fastq-replaceref/$rfq 2> $2/05.fastq-replaceref/$rfq.stderr.log;
		
		mv $2/05.fastq-replaceref/$rfq.clean $2/05.fastq-replaceref/$rfq.tmp
        bbsplitpairs.sh in=$2/05.fastq-replaceref/$rfq.tmp out=$2/05.fastq-replaceref/$rfq.clean outs=$2/05.fastq-replaceref/$rfq.single fint

	fi
	echo -e "\tRead sequences replaced....[x]";


# step-5: Combine, shuffle, separate into R1 and R2 for SRA submission

	f1=$(basename $r1);
	f2=$(basename $r2);

	if [ ! -f $2/06.fastq-forSRA/$f1 ];then
		cat $2/03.fastq-scrubbed/$scrb $2/05.fastq-replaceref/$rfq | sed 's/\sreplaced.*//g' | \
		fasten_randomize -n "${NSLOTS}" -p | \
		fasten_shuffle -n "${NSLOTS}" -d \
		-1 $2/06.fastq-forSRA/$f1 -2 $2/06.fastq-forSRA/$f2;
	fi
	echo -e "\tReads prepared for SRA....[x]";


# step-6: Summarize some count metrics

	sum=$(basename $intlv .fastq).stats.tsv;

	RAW=$(wc -l $2/02.fastq-interleaved/$intlv | awk '{print $1/8}');
        SCRpass=$(wc -l $2/03.fastq-scrubbed/$scrb | awk '{print $1/8}');
        T2Treplaced=$(wc -l $2/05.fastq-replaceref/$rfq | awk '{print $1/8}');
        TOTsra=$(wc -l $2/06.fastq-forSRA/$f1 | awk '{print $1/4}');
        FHS=$(echo "scale=4; $T2Treplaced/$TOTsra" | bc | awk '{printf "%.4f", $0}' );
        SCRfail=$(wc -l $2/03.fastq-scrubbed/$rmvd | awk '{print $1/8}');

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

for i in $(find $1/*R1*fastq); do

	anonymize $i $2 $3;

done

cat $2/07.summary/*tsv | head -1 > $2/Summary.tsv;
cat $2/07.summary/*tsv | grep -v "^ReadID" | sort >> $2/Summary.tsv;

