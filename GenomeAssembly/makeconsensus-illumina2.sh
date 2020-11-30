REF_FILE=$1
BAM_FILE=$2
SEQ_NAME=`basename $REF_FILE _ref.fasta`
OUTPUT_FILE=$3

REF_NAME=`cat $REF_FILE | grep '>' | tr -d '>' | cut -d ' ' -f 1`
LENGTH=`tail -n +2 $REF_FILE | tr -d '\n' | wc -m | xargs`
echo -e "$REF_NAME\t$LENGTH" > my.genome
bedtools bamtobed -i $BAM_FILE > reads.bed
bedtools genomecov -bga -i reads.bed -g my.genome | awk -e '$4 < 1' > zero.bed
maskFastaFromBed -fi $REF_FILE -bed zero.bed -fo masked.fasta
bcftools mpileup -Ou -C50 -f masked.fasta $BAM_FILE | bcftools call --ploidy 1 -mv -Oz -o test.vcf.gz
bcftools index test.vcf.gz
cat masked.fasta | bcftools consensus test.vcf.gz > new_consensus.fasta
echo ">$SEQ_NAME" > $OUTPUT_FILE
tail -n +2 new_consensus.fasta >> $OUTPUT_FILE
