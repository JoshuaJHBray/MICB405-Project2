#!/bin/bash
#prokka_new.sh

#Prepwork for analyses
cd ~
mkdir project2
mkdir project2/prokka_cat
mkdir project2/BWA_output
mkdir project2/RPKM_outputs

#Run Prokka for each MAG .fa

for sample in /projects/micb405/resources/project_2/2018/SaanichInlet_200m/MetaBAT2_SaanichInlet_200m/MedQPlus_MAGs/*.fa
do
	number=$(echo $sample | sed 's/.fa//;s/.*200m.//')
	sid=$(echo $sample | sed 's/.fa//;s/.*\///')
	tax=$(grep -w $sid /projects/micb405/resources/project_2/2018/SaanichInlet_200m/MetaBAT2_SaanichInlet_200m/gtdbtk_output/gtdbtk.*.classification_pplacer.tsv | awk '{ print $2 }' | awk -F";" ' {print $1}' | sed 's/d__//g')
	echo $number, $sid,$tax
	prokka \
	--kingdom $tax \
	--outdir project2/prokka_output/${sid} \
	--force \
	--prefix MAG${number} \
	--cpus 8 \
	$sample	
done

#Concatenate all of the .faa and .fna outputs for KAAS and indexing

cat project2/prokka_output/*/MAG*.faa >> project2/prokka_cat/allMAGs.faa
cat project2/prokka_output/*/MAG*.fna >> project2/prokka_cat/allMAGs.fna

#Index the concatenated reference for the metatranscriptome analysis

bwa index -p ~/project2/BWA_output/allMAGs_index \
/project2/prokka_cat/allMAGs.fna

echo done indexing allMAGs... 

echo Aligning cruise RNA fastqs to the indexed reference

for SI_RNA in /projects/micb405/resources/project_2/2018/Metatranscriptomes/*200m*
do
	cruise=$(echo $SI_RNA | sed 's/_200m.*//;s/.*\///')
	echo Aligning $cruise
	bwa mem -t 8 ~/project2/BWA_output/allMAGs_index \
	$SI_RNA \
	1>~/project2/BWA_output/${cruise}_allMAGs.sam \
	2>~/project2/BWA_output/${cruise}_allMAGs_log.txt 
	echo done aligning ${cruise}
done

#Normalization using RPKM

echo Beginning RPKM analysis...

for SI in ~/project2/BWA_output/*_allMAGs.sam
do
	cruise=$(echo $SI | sed 's/_allMAGs.sam//;s/.*\///')
	echo Running RPKM on $cruise
	/projects/micb405/resources/project_2/2018/rpkm \
	-c ~/project2/prokka_cat/allMAGs.fna \
	-a ~/project2/BWA_output/${cruise}_allMAGs.sam \
	-o ~/project2/RPKM_outputs/${cruise}_allMAGs_RPKM.csv
	echo Done RPKM normalizing ${cruise}
done


