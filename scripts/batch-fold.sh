#!/bin/bash

# Jason Apostol - March 2022
# File to automate the folding of many proteins.
# To be run after download.py, align.py

usage() {
  cat << EOF >&2
Usage: $PROGNAME -f <file> [-o <dir>] [...]
 Required
 -f <file>: the input file containing sequences in FASTA format
 Reccomended
 -a <dir> : the directory for alphafold, default: current dir
 Optional
 -o  <dir>: the output dir, default: data/folded-structures
 -d <dir> : the directory for alphafold's data dir, default: alphafold/data
 -c       : database setting, default: reduced_dbs
 -t       : newest template release date to consider (ISO-8601 format - i.e. YYYY-MM-DD), default: 2022-03-01
 -g       : Use GPU, default true
EOF
  exit 1
}

fasta=""
outdir="../data/folded-structures"
alphadir="."
dbs="reduced_dbs"
t="2022-03-01"
use_gpu=true
gpu_relax=true

while getopts f:o:a:d:c:t:g:r: flag
do
  case $flag in
    (f) fasta=$OPTARG;;
    (o) outdir=$OPTARG ;;
    (a) alphadir=$OPTARG ;;
    (d) datadir=$OPTARG ;;
    (c) dbs=$OPTARG ;;
    (t) t=$OPTARG ;;
    (g) use_gpu=$OPTARG ;;
    (r) gpu_relax=$OPTARG ;;
    (*) usage
       exit 1 ;;
  esac
done

datadir="$alphadir/data"

curr=""
while IFS= read -r line
do
  if [[ ">" == "${line:0:1}" ]]
  then
    curr=${line:1:4}
    echo "----Starting protein $curr----"
  else
    echo ">$curr
$line" > $curr.fasta
    echo $alphadir/run_alphafold.sh -d $datadir -o $outdir -z $alphadir -f ./$curr.fasta -t $t -c $dbs -g $use_gpu -r $gpu_relax
    $alphadir/run_alphafold.sh -d $datadir -o $outdir -z $alphadir -f ./$curr.fasta -t $t -c $dbs -g $use_gpu -r $gpu_relax
    rm $curr.fasta
  fi
done < "$fasta"