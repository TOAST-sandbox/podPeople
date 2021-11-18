#!/usr/bin/env bats

# https://github.com/bats-core/bats-core

load "inc/environment"

function note(){
  echo "# $1" >&3
}

@test "Environment" {
  run bash src/podPeople.sh --check
  echo "$output"
  if [ "$status" -gt 0 ]; then
    echo "ERROR: environment needs work"
  fi
}

