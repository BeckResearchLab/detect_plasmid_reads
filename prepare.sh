#!/bin/bash

PLASMIDS_PATH=/work/data/NCBI_plasmids/plasmid
REFSEQ_PATH=/work/data/refseq
THREADS=24
MAX_PLASMID_LEN=500000
MIN_REFSEQ_LEN=1000000
RANDOM_SEED=42
CLASS_SAMPLES=0

if [ ! -e plasmid_seq.fna ]; then
	./plasmid_seq_extractor.py --plasmids_fna $PLASMIDS_PATH/ncbi_plasmid.fna \
			--plasmids_gbff $PLASMIDS_PATH/ncbi_plasmid.gbff \
			--output_file plasmid_seq.fna --max_length $MAX_PLASMID_LEN
fi

if [ ! -e plasmid_reads1.fq -o ! -e plasmid_reads2.fq ]; then
	# these numbers were taken from an actual sequencing run
	# deviation 105
	# fragment length 625 (150 + 150 + 325)
	/work/software/art/art_illumina --seqSys HS25 --len 150 --paired \
			--sdev 105 --mflen 625 \
			--fcov 1.0 --insRate 0 --insRate2 0 --delRate 0 --delRate2 0 \
			--maxIndel 0 --maskN 0 --noALN --out plasmid_reads \
			--rndSeed $RANDOM_SEED --in plasmid_seq.fna
fi

if [ ! -e plasmid_reads.tsv ]; then
	paste plasmid_reads1.fq plasmid_reads2.fq | awk '{ if (line % 4 == 1) printf("%s%s\t1\n", $1, $2); ++line; }' > plasmid_reads.tsv
fi

if [ ! -e refseq_seq.fna ]; then
	./refseq_seq_extractor.py --refseq_path $REFSEQ_PATH \
			--output_file refseq_seq.fna --threads $THREADS
fi

if [ ! -e refseq_reads1.fq -o ! -e refseq_reads2.fq ]; then
	# these numbers were taken from an actual sequencing run
	# deviation 105
	# fragment length 625 (150 + 150 + 325)
	/work/software/art/art_illumina --seqSys HS25 --len 150 --paired \
			--sdev 105 --mflen 625 \
			--fcov 1.0 --insRate 0 --insRate2 0 --delRate 0 --delRate2 0 \
			--maxIndel 0 --maskN 0 --noALN --out refseq_reads \
			--rndSeed $RANDOM_SEED --in refseq_seq.fna
fi

if [ ! -e refseq_reads.tsv ]; then
	paste refseq_reads1.fq refseq_reads2.fq | awk '{ if (line % 4 == 1) printf("%s%s\t0\n", $1, $2); ++line; }' > refseq_reads.tsv
fi

if [ ! -e all_reads.tsv ]; then
	echo "sequence\tis_plasmid" > all_reads.tsv
	cat plasmid_reads.tsv refseq_reads.tsv >> all_reads.tsv
fi

if [ ! -e balanced_reads.tsv ]; then
	./all_seq_balance.py --input_file all_seq.tsv --output_file balanced_seq.tsv \
			--positive_samples $CLASS_SAMPLES --random_seed $RANDOM_SEED
fi

if [ ! -e all_seq_train.h5 -o ! -e all_seq_valid.h5 -o ! -e all_seq_test.h5 ]; then
	./all_seq_save_hdf5.py --input_file balanced_reads.tsv \
			--train_frac 0.91 --valid_frac 0.01 --test_frac 0.08 \
			--train_file all_seq_train.h5 --valid_file all_seq_valid.h5 \
			--test_file all_seq_test.h5
fi
