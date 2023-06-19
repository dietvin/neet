#!/bin/bash

inbam=/home/vincent/masterthesis/data/curlcake_m6a_rep1.bam
ref=/home/vincent/masterthesis/data/curlcake_ref.fa

#samtools mpileup -f $ref -A -d 0 -Q 0 $inbam > ./curlcake_m6a.pileup

python pileup_extractor.py -i ./curlcake_m6a.pileup -o ./curlcake_m6a_extracted.tsv -r $ref \
	-n 80 \
	-q 10

python neighbour_searcher.py -i ./curlcake_m6a_extracted.tsv -o ./curlcake_m6a_extracted_w_neighbour.tsv -w 2 -t 0.6
