echo """You should modify this section first!"""
files=${1:-"JZLC_04.read1.fq.gz JZLC_04.read2.fq.gz"}
sampleID=${2:-"JZLC_04"}
K=${3:-30}
CPU=${4:-8}
# This depends on your way of installation.
# We created a python virtual environment for metaphlan4;
# modify this if you used a different method.
METAPHLAN=/path/to/metaphlan/bin/activate

# check the UProC installation page to download the database and model
PATH=/path/to/uproc/bin/:$PATH
DB_PFAM=/path/to/pre-compiled/pfam/database
MODEL=/path/to/pre-trained/uproc/model

PATH=/path/to/KMC/bin/:$PATH

echo """===================================="""

echo """annotate using metaphlan4"""
mkdir -p otu
### modify this if you used different method to install metaphlan4
source $METAPHLAN
tmpF=$(echo "$files" | sed 's/ /,/g')
metaphlan --nproc ${CPU} \
  -t rel_ab_w_read_stats \
  --input_type fastq \
  --bowtie2out ${sampleID}.bowtie2.bz2 \
  -o otu/${sampleID}.txt $tmpF
deactivate

echo """annotate using uproc"""
mkdir -p pfam
uproc-dna -t ${CPU} -o pfam/${sampleID}.csv ${DB_PFAM} ${MODEL} ${files}

echo """annotate using kmc"""
mkdir -p K${K}
words=($files)
for word in "${words[@]}";
do
  echo "$word" >> ${sampleID}.lst
done

kmc -k${K} -t${CPU} -cs100000 @${sampleID}.lst $sampleID K${K}
kmc_tools -t${CPU} transform $sampleID dump -s K${K}/${sampleID}.all.txt
# make sure you have enough hard disk to host ${sampleID}.all.txt, as it's going to be huge!
