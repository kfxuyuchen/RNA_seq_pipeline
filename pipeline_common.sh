#!/bin/bash
# Shared helpers for RNA_seq_pipeline (source from other scripts).
# shellcheck shell=bash

# Default PE FASTQ naming (backward compatible with *clean_1.fastq.gz)
: "${R1_SUFFIX:=_1.fastq.gz}"
: "${R2_SUFFIX:=_2.fastq.gz}"
: "${R1_GLOB:=*clean_1.fastq.gz}"

fq_sample_name() {
  # strip directory and R1 suffix → sample id
  local r1="$1"
  local base
  base=$(basename "$r1")
  if [[ "$base" == *"$R1_SUFFIX" ]]; then
    echo "${base%"$R1_SUFFIX"}"
  else
    echo "ERROR: R1 file does not end with R1_SUFFIX='$R1_SUFFIX': $base" >&2
    return 1
  fi
}

fq_mate_path() {
  # given R1 path, return expected R2 path
  local r1="$1"
  local dir base sample
  dir=$(dirname "$r1")
  sample=$(fq_sample_name "$r1") || return 1
  echo "${dir}/${sample}${R2_SUFFIX}"
}

find_r1_files() {
  # print sorted R1 paths under $1
  local dir="$1"
  find "$dir" -maxdepth 1 -name "$R1_GLOB" -print | sort
}

parse_fq_naming_args() {
  # consume one option from "$@"; echo remaining via global __REST
  # usage inside case: --r1-suffix) R1_SUFFIX="$2"; shift 2 ;;
  :
}
