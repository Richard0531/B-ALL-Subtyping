r1=$1
r2=$2
samplename=$3
outfolder=$4
rundate=$5
#module load star/2.7.9a
source "/home/hcliulab/miniconda3/etc/profile.d/conda.sh"
conda activate subtyping
genomeDir=/home/hcliulab/yhsiao/GRCh38/STAR_index/
hg38_ref=/home/hcliulab/yhsiao/GRCh38/GRCh38.primary_assembly.genome.fa
PICARD=/home/hcliulab/yhsiao/GRCh38/picard.jar
tmpdir=/home/hcliulab/yhsiao/GRCh38/tmp/
cd $outfolder
mkdir $samplename
cd $samplename
mkdir 1pass
cd 1pass

readCommandParam=""
if [[ "$r1" == *.gz ]]
then
	readCommandParam="--readFilesCommand zcat"
fi

STAR $readCommandParam --genomeDir $genomeDir --readFilesIn $r1 $r2 --runThreadN 24 --chimSegmentMin 12 --chimJunctionOverhangMin 12 --chimOutType Junctions --chimOutJunctionFormat 1
genomeDir2=$outfolder/$samplename/ref_2pass
cd ../
mkdir ref_2pass
cd ref_2pass
STAR --runMode genomeGenerate --genomeDir $genomeDir2 --genomeFastaFiles $hg38_ref \
     --sjdbFileChrStartEnd $outfolder/$samplename/1pass/SJ.out.tab --sjdbOverhang 75 --runThreadN 24
cd ../
mkdir 2pass
cd 2pass
STAR $readCommandParam --genomeDir $genomeDir2 --readFilesIn $r1 $r2 --runThreadN 24 --chimSegmentMin 12 --chimJunctionOverhangMin 12 --chimOutType Junctions --chimOutJunctionFormat 1

conda deactivate
### add read group and sort
java -Djava.io.tmpdir=$tmpdir -Xmx40g -jar $PICARD AddOrReplaceReadGroups -I Aligned.out.sam -O Aligned.rg_added_sorted.bam -SO coordinate -RGID ${samplename} -RGLB ${samplename} -RGPL Illumina -RGPU Next_Hi_Seq  -RGSM ${samplename} -DT ${rundate}

### markduplicate
java -Djava.io.tmpdir=$tmpdir -Xmx40g -jar $PICARD MarkDuplicates -I Aligned.rg_added_sorted.bam  -O ${samplename}.rg_added_sorted.dupmarked.bam  -CREATE_INDEX true -VALIDATION_STRINGENCY SILENT -M output.metrics

rm Aligned.out.sam
rm Aligned.rg_added_sorted.bam
cd ../
#rm -R 1pass
#rm -R ref_2pass
mv 2pass/${samplename}.rg_added_sorted.dupmarked.ba* .
mv 2pass/Chimeric.out.junction Chimeric.out.junction
mv 2pass/output.metrics output.metrics
mv 2pass/SJ.out.tab SJ.out.tab
mv 2pass/Log.final.out Log.final.out
#rm -R 2pass

conda activate subtyping
folder=`dirname $input`
featureCounts -s 2 -p -B -a /home/hcliulab/yhsiao/GRCh38/STAR_index/gencode.v39.primary_assembly.annotation.gtf -o $outfolder/$samplename/count_unique.txt  ${samplename}.rg_added_sorted.dupmarked.bam
featureCounts -O -M -p -B -T 10 -a /home/hcliulab/yhsiao/GRCh38/STAR_index/gencode.v39.primary_assembly.annotation.gtf -o $outfolder/$samplename/${samplename}_count_multi.txt ${samplename}.rg_added_sorted.dupmarked.bam
cat $outfolder/$samplename/${samplename}_count_multi.txt | tail -n+3 | cut -f1,7 | sed -e "s|\.[0-9]*||"  > $outfolder/$samplename/${samplename}_count_multi_simplified.txt
cat $outfolder/$samplename/${samplename}_count_multi.txt | tail -n+3 | cut -f1,7   > $outfolder/$samplename/${samplename}_count_multi_simplified_2.txt

Rscript /home/hcliulab/yhsiao/gene_expression/RNA_pipeline_Liu_MDALL.R $samplename $outfolder
#Rscript /home/hcliulab/yhsiao/gene_expression/RNA_pipeline_Liu_T.R $samplename $outfolder
STAR-Fusion --genome_lib_dir /home/hcliulab/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ -J Chimeric.out.junction --output_dir $outfolder/${samplename}

source "$HOME/.sdkman/bin/sdkman-init.sh"
mkdir $outfolder/$samplename/variant_calling

base=$samplename.rg_added_sorted.dupmarked
GATK=/home/hcliulab/yhsiao/GATK/GenomeAnalysisTK.jar
hg38_ref=/home/hcliulab/yhsiao/GRCh38/GRCh38.primary_assembly.genome.fa
java -Djava.io.tmpdir=/home/hcliulab/yhsiao/tmp/  -Xmx40g -jar $GATK -T SplitNCigarReads -R $hg38_ref -I ${base}.bam -o $outfolder/$samplename/variant_calling/${base}.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
cd variant_calling
cat /home/hcliulab/yhsiao/variant_calling/chromosomes | xargs -I {} -P 24 java -Djava.io.tmpdir=/home/hcliulab/yhsiao/tmp/ -Xmx8g -jar $GATK  -T HaplotypeCaller -L {} -R $hg38_ref -I ${base}.split.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o raw.scatter{}.vcf
java -cp $GATK org.broadinstitute.gatk.tools.CatVariants \
	-R $hg38_ref \
        -V raw.scatterchr1.vcf \
        -V raw.scatterchr2.vcf \
        -V raw.scatterchr3.vcf \
        -V raw.scatterchr4.vcf \
        -V raw.scatterchr5.vcf \
        -V raw.scatterchr6.vcf \
        -V raw.scatterchr7.vcf \
        -V raw.scatterchr8.vcf \
        -V raw.scatterchr9.vcf \
        -V raw.scatterchr10.vcf \
        -V raw.scatterchr11.vcf \
        -V raw.scatterchr12.vcf \
        -V raw.scatterchr13.vcf \
        -V raw.scatterchr14.vcf \
        -V raw.scatterchr15.vcf \
        -V raw.scatterchr16.vcf \
        -V raw.scatterchr17.vcf \
        -V raw.scatterchr18.vcf \
        -V raw.scatterchr19.vcf \
        -V raw.scatterchr20.vcf \
        -V raw.scatterchr21.vcf \
        -V raw.scatterchr22.vcf \
        -V raw.scatterchrX.vcf \
        -V raw.scatterchrY.vcf \
        -out ${base}.vcf \
        -assumeSorted
java -Djava.io.tmpdir=/home/hcliulab/yhsiao/tmp/  -jar $GATK -T VariantFiltration -R $hg38_ref -V ${base}.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ${base}.filtered.vcf
cat ${base}.filtered.vcf | grep "PASS" > ${base}.filtered.PASSonly.vcf
vep --fork 12 --offline -i ${base}.filtered.vcf -o ${base}.filtered.vep.vcf  --dir_cache /home/hcliulab/yhsiao/vc/ --polyphen b --sift b --symbol --cache_version 113 --check_existing --af --max_af --af_1kg --ccds --canonical --protein --vcf
filter_vep --input_file ${base}.filtered.vep.vcf  --filter "(CANONICAL is YES) and IMPACT in HIGH,MODERATE" --only_matched > ${base}.filtered.PASSonly.vep.coding.vcf
gzip ${base}.filtered.PASSonly.vcf

cd ..
mkdir $outfolder/$samplename/cnv
cd cnv
vcf=$outfolder/$samplename/variant_calling/${base}.filtered.vcf
count_file=$outfolder/$samplename/${samplename}_count_multi_simplified.txt
output_folder=$outfolder/$samplename/cnv

random_id=`echo $RANDOM | md5sum | head -c20`
vcf_folder=`dirname $vcf`
vcf_name=`basename $vcf`
cat $count_file > ${random_id}_count_file
echo "out_dir = \"$output_folder\"" > ${random_id}.config
echo "out_dir = \"$output_folder\"" > ${random_id}.config
echo "count_dir = \"$output_folder\"" >> ${random_id}.config
echo "snv_dir = \"${vcf_folder}/\"" >> ${random_id}.config
echo -e "${samplename}\t${random_id}_count_file\t$vcf_name" > ${random_id}.meta
Rscript -e "library(RNAseqCNV); library(showtext); showtext_auto(); RNAseqCNV_wrapper(config = \"${random_id}.config\", metadata = \"${random_id}.meta\", snv_format = 'vcf')"
rm ${random_id}_count_file
rm ${random_id}.config
rm ${random_id}.meta










