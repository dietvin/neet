#!/bin/bash

ref=/home/vincent/masterthesis/data/curlcake_ref.fa


inbam=/home/vincent/masterthesis/data/curlcake_m6a_rep1.bam
inbam=/home/vincent/masterthesis/data/curlcake_unm_rep1.bam

samtools mpileup -f $ref -A -d 0 -Q 0 $inbam > ./curlcake_m6a.pileup
samtools mpileup -f $ref -A -d 0 -Q 0 $inbam > ./curlcake_unm.pileup

python pileup_extractor.py -i ./curlcake_m6a.pileup -o ./curlcake_m6a_extracted.tsv -r $ref \
 	-n 80 \
 	-q 10

python pileup_extractor.py -i ./curlcake_unm.pileup -o ./curlcake_unm_extracted.tsv -r $ref \
	-n 80 \
	-q 10

python neighbour_searcher.py -i ./curlcake_m6a_extracted.tsv -o ./curlcake_m6a_extracted_w_neighbour.tsv -w 2 -t 0.6

python neighbour_searcher.py -i ./curlcake_unm_extracted.tsv -o ./curlcake_unm_extracted_w_neighbour.tsv -w 2 -t 0.6






