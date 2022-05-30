mkdir ./tempdir
export TMPDIR=./tempdir
inbam=$1 ### bam file
intsv=$2 ### cell bin segmentation file with coordinates and cell ID
outbamdir=$3 ### out bam path
perfix=$4   ### sample ID
outloomdir=$5 ### output loom path

perl change_CB.pl -i $intsv -b $inbam > ${outbamdir}/${perfix}.sam
samtools view -Sb ${outbamdir}/${perfix}.sam -o ${outbamdir}/${perfix}.bam
samtools sort ${outbamdir}/${perfix}.bam -o ${outbamdir}/${perfix}.sort.bam -@ 3 -m 60000M
samtools index ${outbamdir}/${perfix}.sort.bam -@ 3
samtools sort -t CB -O BAM -o ${outbamdir}/cellsorted_${perfix}.sort.bam ${outbamdir}/${perfix}.sort.bam -@ 3 -m 60000M
velocyto run ${outbamdir}/${perfix}.sort.bam AmexT_v47.ALL.FINAL.annotated_0925_renameID.gtf --outputfolder $outloomdir --samtools-threads 20 --samtools-memory 60000
