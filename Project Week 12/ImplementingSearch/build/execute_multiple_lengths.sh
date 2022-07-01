#!/bin/bash
for len in 40 60 80 100 
do 
  echo "Benchmarking for Length $len"
  /usr/bin/time -v ./bin/fmindex_search --index myIndex.index --query ../data/illumina_reads_$len.fasta.gz &> fmindex_q100000_l$len.log
  egrep 'The Program took|bytes' fmindex_q100000_l$len.log 
done 