#!/usr/bin/env bash


source /opt/sge/default/common/settings.sh
source /etc/profile.d/modules.sh
source ~/.bashrc

module load \
 BBMap/38.90 \
 BEDTools/2.27.1 \
 bowtie2/2.3.5.1 \
 gcc/4.9.3 \
 Mash/2.0 \
 samtools/1.10

export LD_LIBRARY_PATH=/apps/x86_64/gcc/4.9.3/lib64:/apps/x86_64/gcc/4.9.3/lib:$LD_LIBRARY_PATH

cd "${o}"


	readid=$(basename $r1 | sed 's/_L001.*//g');

# step-1: Interleave R1 and R2

	r2=$(echo $r1 | sed 's/_R1_/_R2_/g');
	intlv=$(basename $r1 | sed 's/R1_001/interleaved/g');

	echo -e "$r1\t$r2\t$intlv";

	if [ ! -f ./02.fastq-interleaved/$intlv ];then
		/scicomp/groups-pure/Projects/TOAST/fasten/target/release/fasten_shuffle -n "${NSLOTS}" -1 $r1 -2 $r2 > \
		./02.fastq-interleaved/$intlv;
	fi

# step-2: Human read scrubber and repair broken pairs

	scrb=$(basename $intlv .fastq)_scrubbed.fastq;
	tmp=$(echo $scrb | sed 's/scrubbed/tmp/g');
	single=$(echo $scrb |sed 's/scrubbed/single/g');
        rtmp=$(echo $scrb | sed 's/scrubbed/rtmp/g');
        rsingle=$(echo $scrb |sed 's/scrubbed/rsingle/g');
        rmvd=$(echo $scrb |sed 's/scrubbed/removed/g');
	echo $scrb;

	if [ ! -f ./03.fastq-scrubbed/$scrb ];then
		/scicomp/groups-pure/Projects/TOAST/sra-human-scrubber/scripts/scrub.sh -r ./02.fastq-interleaved/$intlv;
		cat ./02.fastq-interleaved/$intlv.clean | sed 's/\s/\//g' | sed 's/\:N\:0.*//g' > ./03.fastq-scrubbed/$tmp;
		cat ./02.fastq-interleaved/$intlv.removed_spots | sed 's/\s/\//g' | sed 's/\:N\:0.*//g' > ./03.fastq-scrubbed/$rtmp;

		bbsplitpairs.sh in=03.fastq-scrubbed/$tmp out=03.fastq-scrubbed/$scrb outs=03.fastq-scrubbed/$single fint
                bbsplitpairs.sh in=03.fastq-scrubbed/$rtmp out=03.fastq-scrubbed/$rmvd outs=03.fastq-scrubbed/$rsingle fint

	fi

# step-3: Bowtie2 map to human reference
		
	bm=$(basename $intlv .fastq).t2t.bam;
	echo $bm;

	if [ ! -f ./04.bam-t2t/$bm ];then
		bowtie2 \
		-X 1000 \
		-p "${NSLOTS}" -x /scicomp/reference/bowtie2/t2t-chm13/t2t-chm13.fasta \
                --interleaved ./03.fastq-scrubbed/$rmvd | samtools view -b -F 12 > ./04.bam-t2t/$bm;

	fi

# step-4: Replace patient sequence with corresponding T2T reference sequence

	rfq=$(basename $bm bam)fastq;
	echo $rfq;

	if [ ! -f ./05.fastq-replaceref/$rfq ];then
		/scicomp/groups-pure/Projects/TOAST/MWeigand/lskScripts/scripts/replaceReadsWithReference.pl \
		/scicomp/reference/bowtie2/t2t-chm13/t2t-chm13.fasta ./04.bam-t2t/$bm 1> ./05.fastq-replaceref/$rfq 2> ./05.fastq-replaceref/$rfq.stderr.log;
		
		mv ./05.fastq-replaceref/$rfq.clean ./05.fastq-replaceref/$rfq.tmp
                bbsplitpairs.sh in=./05.fastq-replaceref/$rfq.tmp out=./05.fastq-replaceref/$rfq.clean outs=./05.fastq-replaceref/$rfq.single fint

	fi

# step-5: Combine, shuffle, separate into R1 and R2 for SRA submission

	f1=$(basename $r1);
	f2=$(basename $r2);

	if [ ! -f ./06.fastq-forSRA/$f1 ];then
		cat ./03.fastq-scrubbed/$scrb ./05.fastq-replaceref/$rfq | sed 's/\sreplaced.*//g' | \
		/scicomp/groups-pure/Projects/TOAST/fasten/target/release/fasten_randomize -n "${NSLOTS}" -p | \
		/scicomp/groups-pure/Projects/TOAST/fasten/target/release/fasten_shuffle -n "${NSLOTS}" -d \
		-1 ./06.fastq-forSRA/$f1 -2 ./06.fastq-forSRA/$f2;
	fi


# step-6: Summarize some count metrics

	sum=$(basename $intlv .fastq).stats.tsv;

	RAW=$(wc -l ./02.fastq-interleaved/$intlv | awk '{print $1/8}');
        SCRpass=$(wc -l ./03.fastq-scrubbed/$scrb | awk '{print $1/8}');
        T2Treplaced=$(wc -l ./05.fastq-replaceref/$rfq | awk '{print $1/8}');
        TOTsra=$(wc -l ./06.fastq-forSRA/$f1 | awk '{print $1/4}');
        FHS=$(echo "scale=4; $T2Treplaced/$TOTsra" | bc | awk '{printf "%.4f", $0}' );
        SCRfail=$(wc -l ./03.fastq-scrubbed/$rmvd | awk '{print $1/8}');
        SCRfinal=$(wc -l ./10.rescrub/$mxscrb | awk '{print $1/8}');
        SCRreplaced=$(wc -l ./05.fastq-replaceref/$rfq.clean | awk '{print $1/8}');

        if [ "$RAW" -gt "$SCRpass" ];then
                T2F=$(echo "scale=4; $T2Treplaced/$SCRfail" | bc | awk '{printf "%.4f", $0}' );
        else
                T2F=0;
        fi

        echo -e "ReadID\tRaw_pairs\tScrub_pass\tScrub_remove\tReplaced_pairs\tReplaced_frac\tTotal_pairs_final\tHuman_frac_final\tScrub_pass_replaced\tScrub_pass_final" > ./07.summary/$sum;
        echo -e "$readid\t$RAW\t$SCRpass\t$SCRfail\t$T2Treplaced\t$T2F\t$TOTsra\t$FHS\t$SCRreplaced\t$SCRfinal" >> ./07.summary/$sum;


echo "INFO: Step 6 completed. Pipeline finished"

