#!/bin/bash

# Jason Apostol - March 2022
# File to automate the folding of many proteins
# To be run after download.py, align.py

usage() {
  cat << EOF >&2
Usage: $PROGNAME -f <file> [-o <dir>] [...]

 -f <file>: the input file containing sequences in FASTA format
 -o  <dir>: the output dir, default: data/folded-structures
 -a <dir> : the directory for alphafold, default: current dir
 -d <dir> : the directory for alphafold's data dir, default: alphafold/data
 -c       : database setting, default: reduced_dbs
 -t       : newest template release date to consider (ISO-8601 format - i.e. YYYY-MM-DD), default: 2022-03-01
EOF
  exit 1
}

fasta=""
outdir="."
alphadir=""
datadir="$alphadir/data"
dbs="reduced_dbs"
t="2022-03-01"

while getopts f:o:a:d:c:t: flag
do
  case $flag in
    (f) fasta=$OPTARG;;
    (o) outdir=$OPTARG ;;
    (a) alphadir=$OPTARG ;;
    (d) datadir=$OPTARG ;;
    (c) dbs=$OPTARG ;;
    (t) t=$OPTARG ;;
    (*) usage
       exit 1 ;;
  esac
done

curr=""
while IFS= read -r t
do
  if [[ ">" == "${t:0:1}" ]]
  then
    curr=${t:1:4}
    echo "----Starting protein $curr----"
  else
    echo ">$curr\n$t" > tmp.fasta
    source "$alphadir/run_alphafold.sh" -d $datadir -o $outdir -f tmp.fasta -t $t -c $dbs
    rm tmp.fasta
  fi
done < "$fasta"