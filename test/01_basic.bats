#!/usr/bin/env bats

# https://github.com/bats-core/bats-core

load "inc/environment"

function note(){
  echo "# $1" >&3
}

@test "Environment" {
  run bash src/podPeople.sh --check
  note "$output"
  if [ "$status" -gt 0 ]; then
    note "ERROR: environment needs work"
    exit $status
  fi
}

@test "basic data" {
  run bash src/podPeople.sh test basic.out /scicomp/reference/bowtie2/t2t-chm13/t2t-chm13.fasta
  note "$output"
  if [ $status -gt 0 ]; then
    note "ERROR with basic test"
    exit 1
  fi
  ls -R basic.out
  column -ts $'\t' basic.out/Summary.tsv
}

