#/bin/bash

REF=ref/chrM-rsrs.fa
REF=ref/rsrs.fasta
REF_TMP=ref/chrM.fa
INFILE=$1
OUTFILE=$1.bam

echo Indexing ref...
cp $REF $REF_TMP
time bwa index $REF_TMP
echo Aligning...
#NOTE: this aligns near 310 in an unpredictable way
time bwa mem -B 4 -t 6 $REF_TMP $INFILE| samtools view -@ 6 -b - > $OUTFILE.unsorted
echo Sorting...
time samtools sort -@ 6 $OUTFILE.unsorted -o $OUTFILE
rm $OUTFILE.unsorted

samtools index $OUTFILE
./bam_to_snp.py $OUTFILE -i 38 --no-y --no-ystr
