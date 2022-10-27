#
# Copyright 2021 Simone Maestri. All rights reserved.
# Simone Maestri <simone.maestri@univr.it>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

#!/bin/bash


PIPELINE_DIR=$(realpath $( dirname "${BASH_SOURCE[0]}" ))
source $PIPELINE_DIR"/config_Variant_calling_pipeline.sh"

usage="$(basename "$0") [-1 reads_R1] [-2 reads_R2] [-r reference]"

while :
do
    case "$1" in
      -h | --help)
          echo $usage
          exit 0
          ;;
      -1)
          reads_R1=$(realpath $2)
          shift 2
          echo "Reads R1: $reads_R1"
          ;;

      -2)
           reads_R2=$2
           shift 2
           echo "Reads_R2: $reads_R2"
           ;;
      -r)
          reference=$2
          shift 2
          echo "Reference: $reference"
           ;;
       --) # End of all options
           shift
           break
           ;;
       -*)
           echo "Error: Unknown option: $1" >&2
           ## or call function display_help
           exit 1
           ;;
        *) # No more options
           break
           ;;
    esac
done

WORKING_DIR=$(realpath $(dirname $reads_R1))
QC=$WORKING_DIR"/QC"
SAMPLE_NAME=$(echo $(basename $reads_R1) | sed 's/_R1.*//')
REF_NAME=$(echo $(basename $reference) | sed 's/\.fa.*//')
BAM=$WORKING_DIR"/"$SAMPLE_NAME"_mapped_to_"$REF_NAME"_tmp1.bam"
BAM_MD=$WORKING_DIR"/"$SAMPLE_NAME"_mapped_to_"$REF_NAME"_tmp2.bam"
BAM_MD_CLIPPED=$WORKING_DIR"/"$SAMPLE_NAME"_mapped_to_"$REF_NAME"_MarkDup_Clipped.bam"

#Reads quality check
$FASTQC $reads_R1 $reads_R2

#Reads trimming
$FASTP -i $reads_R1 -I $reads_R2 -o $SAMPLE_NAME"_trimmed1.fastq.gz" -O $SAMPLE_NAME"_trimmed2.fastq.gz" --failed_out $SAMPLE_NAME"_failed.fastq.gz" \
--trim_poly_x  --overrepresentation_analysis -h $SAMPLE_NAME"_adapters_removal_report.html" -w $THREADS

#Indexing reference
if [ ! -f "$reference.sa" ]; then
  $BWA index $reference
fi

if [ ! -f "$reference.fai" ]; then
  $SAMTOOLS faidx $reference
fi

if [ ! -f "$reference.dict" ]; then
  $PICARD CreateSequenceDictionary R=$reference O=$reference".dict"
fi

#Perform alignment
$BWA mem -R $(echo "@RG\tID:$SAMPLE_NAME\tPU:lane\tLB:$SAMPLE_NAME\tSM:$SAMPLE_NAME\tCN:Illumina_Variant_calling\tPL:ILLUMINA") -t $THREADS $reference $SAMPLE_NAME"_trimmed1.fastq.gz" $SAMPLE_NAME"_trimmed2.fastq.gz" | $SAMTOOLS sort --threads $THREADS -o $BAM

#Index BAM
$SAMTOOLS index $BAM

#Statistics
$SAMTOOLS flagstat $BAM > $BAM.flagstat
$SAMTOOLS stats $BAM > $BAM.stats
$PICARD CollectInsertSizeMetrics I=$BAM O=$SAMPLE_NAME"_tmp1_insert_size_metrics.txt" H=$SAMPLE_NAME"_tmp1_insert_size_histogram.pdf"

#MarkDuplicates
$PICARD MarkDuplicates I=$BAM O=$BAM_MD M=$SAMPLE_NAME"_MarkDuplicates_metrics.txt" REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=false
$SAMTOOLS index $BAM_MD

#Delete tmp BAM file
rm $BAM $BAM.bai

#ClipOverlaps
$FGBIO ClipBam -i $BAM_MD -o $BAM_MD_CLIPPED -r $reference --clip-overlapping-reads -c Hard -m $SAMPLE_NAME"_fgbio_metrics.txt"

#Delete tmp BAM file
rm $BAM_MD $BAM_MD.bai

#Statistics
$SAMTOOLS flagstat $BAM_MD_CLIPPED > $BAM_MD_CLIPPED.flagstat
$SAMTOOLS stats $BAM_MD_CLIPPED > $BAM_MD_CLIPPED.stats
$PICARD CollectInsertSizeMetrics I=$BAM_MD_CLIPPED O=$SAMPLE_NAME"_insert_size_metrics.txt" H=$SAMPLE_NAME"_insert_size_histogram.pdf"

#Perform variant calling with HaplotypeCaller
$GATK HaplotypeCaller -R $reference -I $BAM_MD_CLIPPED -ERC GVCF --output $SAMPLE_NAME".complete.raw.g.vcf" --standard-min-confidence-threshold-for-calling 30.0 --dont-use-soft-clipped-bases true --sample-ploidy $PLOIDY

#Genotype gVCF
$GATK GenotypeGVCFs -R $reference -V $SAMPLE_NAME".complete.raw.g.vcf" -G StandardAnnotation -O $SAMPLE_NAME".complete.raw.vcf"

#Index featurefile
$GATK IndexFeatureFile -I $SAMPLE_NAME".complete.raw.vcf"

#Filter SNPs
$GATK SelectVariants --select-type-to-include SNP --output $SAMPLE_NAME".snps.raw.vcf" -V $SAMPLE_NAME".complete.raw.vcf"

$GATK VariantFiltration -R $reference -V $SAMPLE_NAME".snps.raw.vcf" \
--filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filter-name "Broad_SNP_filter" -O $SAMPLE_NAME".snps.filtered.vcf" ;

#Delete unfiltered SNPs file
rm $SAMPLE_NAME".snps.raw.vcf" rm $SAMPLE_NAME".snps.raw.vcf.idx"

#Filter InDels
$GATK SelectVariants --select-type-to-exclude SNP --output $SAMPLE_NAME".indels.raw.vcf" -V $SAMPLE_NAME".complete.raw.vcf"

$GATK VariantFiltration -R $reference -V $SAMPLE_NAME".indels.raw.vcf" \
--filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
--filter-name "Broad_indel_Filter" -O $SAMPLE_NAME".indels.filtered.vcf" ;

#Delete unfiltered InDels file
rm $SAMPLE_NAME".indels.raw.vcf" $SAMPLE_NAME".indels.raw.vcf.idx"

#Remove complete unfiltered vcf
rm $SAMPLE_NAME".complete.raw.vcf" $SAMPLE_NAME".complete.raw.vcf.idx"

#Merge SNPs and InDels
$GATK MergeVcfs -I  $SAMPLE_NAME".snps.filtered.vcf" -I $SAMPLE_NAME".indels.filtered.vcf" -O $SAMPLE_NAME".variants.filtered_tmp.vcf"
bgzip $SAMPLE_NAME".variants.filtered_tmp.vcf"
tabix $SAMPLE_NAME".variants.filtered_tmp.vcf.gz"

#Keep only filtered variants
$GATK SelectVariants -R $reference --variant $SAMPLE_NAME".variants.filtered_tmp.vcf.gz" --exclude-filtered -O $SAMPLE_NAME".variants.filtered.vcf"
bgzip $SAMPLE_NAME".variants.filtered.vcf"
tabix $SAMPLE_NAME".variants.filtered.vcf.gz"

#Compress gVCF
bgzip $SAMPLE_NAME".complete.raw.g.vcf"
tabix $SAMPLE_NAME".complete.raw.g.vcf.gz"

#Delete tmp files
rm $SAMPLE_NAME".variants.filtered_tmp.vcf.gz" $SAMPLE_NAME".variants.filtered_tmp.vcf.gz.tbi" $SAMPLE_NAME".variants.filtered_tmp.vcf.idx" $SAMPLE_NAME".indels.filtered.vcf" $SAMPLE_NAME".snps.filtered.vcf" \
$SAMPLE_NAME".indels.filtered.vcf.idx" $SAMPLE_NAME".snps.filtered.vcf.idx"  $SAMPLE_NAME".variants.filtered.vcf.idx" $SAMPLE_NAME".complete.raw.g.vcf.idx"

#Run Qualimap
$QUALIMAP bamqc -bam $BAM_MD_CLIPPED -c -nt $THREADS --java-mem-size=500G

#Move quality check files to QC directory
mkdir $QC
mv *html *stat* *pdf *txt *zip *json $QC
