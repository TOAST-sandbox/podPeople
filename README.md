# Read_anonymizer

[[_TOC_]]

## Overview  
This simple workflow to replaces the sequence of human read contaminants from individuals (patients) with the corresponding sequence in a published human reference genome. Motivation for this project came from a need to prepare technical benchmark datasets containing human contamination with:  

1. known quantities,
1. more realistic distribution than simulated reads, 
1. read sequences derived from published references

---

```mermaid

graph TD
  R1(R1.fastq):::classfi --- SH[[fasten_shuffle]] --> IN(interleaved.fastq)
  R2(R2.fastq):::classfi --> SH
  IN --> SCRUB[[SRA Scrubber]]
  SCRUB --- BB1[[bbsplitpairs.sh]] --> CLEAN(clean.fastq)
  SCRUB --- BB2[[bbsplitpairs.sh]] --> REM(removed.fastq)
  REM --> BOW[[bowtie2]]
  HS(Human ref genome):::classfi --> BOW
  BOW --> BAM(mapped.bam) 
  BAM --> REP[[replaceReadsWithReference.pl]]
  HS --> REP
  REP --- BB3[[bbsplitpairs.sh]]--> READS(replaced.fastq)
  READS --> COMB[[fasten_randomize]]
  CLEAN --> COMB
  COMB --- SHUF[[fasten_shuffle]] --> F1(final.R1.fastq):::classfo
  SHUF --> F2(final.R2.fastq):::classfo

  classDef classfi fill:#85C1E9
  classDef classfo fill:#52BE80

```
---
## Workflow scripts

### Internal use only

1. `ReadCleanBot.internal.sh`  
Pipeline run within TOAST private scicomp workspace to prepare reads for dataset 6 of Xiaoli *et al.*.

1. `ReadCleanBot.aspen.sh`  
Executable for submitting `ReadCleanBot.internal.sh` to Aspen HPC for improved runtime.  

### Generalized for distribution

1. `ReadCleanBot.general.sh`  
Simplified and generalized for external distribution.  
 
    *Command:* 
    ~~~
    ./src/ReadCleanBot.general [in-dir] [out-dir] [path-to-ref]
        [in-dir]      Path to directory of input reads pairs in fastq format with filenames like `name_S1_L001_R[12]_001.fastq`
        [out-dir]     Path to new output directory.
        [path-to-ref] Path to human reference genome assembly fasta (must be indexed for bowtie2).
    ~~~


## Dependencies
This workflow uses:
1. NCBI human read removal tool (aka 'SRA scrubber') -- https://github.com/ncbi/sra-human-scrubber  
1. Bowtie2  
1. ReplaceReadsWithReference.pl -- https://github.com/lskatz/lskScripts/blob/master/scripts/replaceReadsWithReference.pl  
1. Fasten -- https://github.com/lskatz/fasten  
1. T2T human reference genome assembly -- https://www.ncbi.nlm.nih.gov/bioproject/559484  
1. BBTools (bbsplitpairs.sh) -- https://jgi.doe.gov/data-and-tools/bbtools/  
