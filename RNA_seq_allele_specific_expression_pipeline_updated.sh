#new version of the pipeline as of 25/03/18 

#I need to pick wither gtf or gff! 

name=iyBomPrat1.1
ref_GCA=GCA_930367275.1
alt_GCA=GCA_930367225.1
fastq_ID=ERR7113577

#RENAME FASTA FILES SO IT HAS CONSISTENT NAMING WITH THE ANNOTATION FILE (GTF)
ref_fasta=/Volumes/Seagate/Frankie_DTOL_lep_project/phase_two_fasta/${name}_${ref_GCA}.fa
ref_fasta_num=/Volumes/Seagate/Frankie_DTOL_lep_project/phase_two_fasta/${name}_${ref_GCA}_num.fa
alt_fasta=/Volumes/Seagate/Frankie_DTOL_lep_project/phase_two_fasta/${name}_${alt_GCA}.fa
alt_fasta_num=/Volumes/Seagate/Frankie_DTOL_lep_project/phase_two_fasta/${name}_${alt_GCA}_num.fa

gunzip $ref_fasta.gz && python /Users/frankieswift/Desktop/Phase_2/name_converter_fasta.py -fa $ref_fasta -out_fa $ref_fasta_num -contig_swap && samtools index $ref_fasta_num && gzip $ref_fasta

gunzip $alt_fasta.gz && python /Users/frankieswift/Desktop/Phase_2/name_converter_fasta.py -fa $alt_fasta -out_fa $alt_fasta_num -long_parse && samtools index $alt_fasta_num && gzip $alt_fasta

#CREATE AN ALTERNATE VCF FILE IF DOESN'T EXIST
alt_fasta_num=/Volumes/Seagate/Frankie_DTOL_lep_project/phase_two_fasta/${name}_${alt_GCA}_num.fa
ref_output=/Volumes/Seagate/Frankie_DTOL_lep_project/phase_two_fasta/${name}_${ref_GCA}.fa
output=/Volumes/Seagate/Frankie_DTOL_lep_project/outputs/bam/${name}_alt_ref_alignment.bam
samtools_output=/Volumes/Seagate/Frankie_DTOL_lep_project/outputs/samtools/${name}_alt_ref_alignment.sort.bam
vcf_output=/Volumes/Seagate/Frankie_DTOL_lep_project/outputs/VCF/${name}_alt_ref.vcf.gz

#all the commands together!
minimap2 -ax asm10 ${alt_fasta_num} ${ref_output} | samtools view -Sb > ${output} && samtools sort ${output} > ${samtools_output} && samtools index ${samtools_output} && bcftools mpileup -f ${alt_fasta_num} --min-MQ 60 --annotate FORMAT/DP ${samtools_output} | bcftools call -m --ploidy 1 | bgzip > ${vcf_output} && bcftools index -f ${vcf_output}

#DOWNLOAD THE PAIRED END RNA SEQ READS - !NEED TO ADD THE LOCATION TO DOWNLOAD THEM TO!
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR628/006/ERR6286706/ERR6286706_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR628/006/ERR6286706/ERR6286706_1.fastq.gz


#ASSESS THE QUAILITY OF THE READS
fastq1=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/${fastq_ID}_1.fastq
fastq2=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/${fastq_ID}_2.fastq

fastqc --nogroup $fastq1 $fastq2

#DO BASE TRIMMING AND LOW QUALITY (PHRED SCORE) BASE REMOVAL
fastq1_base_qual_check=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/fastp_output/${fastq_ID}_1_base_quaility_checked.fastq
fastq2_base_qual_check=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/fastp_output/${fastq_ID}_2_base_quaility_checked.fastq
fastp_html_report=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/fastp_report/${name}_${fastq_ID}_fastp_report.html
fastp_json_report=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/fastp_report/${name}_${fastq_ID}_fastp_report.json

fastp -i $fastq1 -I $fastq2 -o $fastq1_base_qual_check -O $fastq2_base_qual_check -q 20 -l 30 -h $fastp_html_report -j $fastp_json_report

#REASSESS THE QUAILITY OF THE READS AND CHECK THEY'RE LOOKING OKAY!
fastqc --nogroup $fastq1_base_qual_check $fastq2_base_qual_check

#CREATE STAR REFERENCE GENOMES FOR BOTH THE REFERENCE AND ALTERNATE ASSEMBILIES
ref_fasta_alt=/Volumes/Seagate/Frankie_DTOL_lep_project/phase_two_fasta/${name}_${ref_GCA}_num.fa
alt_fasta_alt=/Volumes/Seagate/Frankie_DTOL_lep_project/phase_two_fasta/${name}_${alt_GCA}_num.fa
ref_gtf=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/Bombus_pratorum-${ref_GCA}-2022_03-genes.gtf
alt_gtf=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/Bombus_pratorum-${alt_GCA}-2022_03-genes.gtf

mkdir -p /Volumes/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/star_index_${name}_ref/
mkdir -p /Volumes/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/star_index_${name}_alt/

STAR --runMode genomeGenerate --genomeDir /Volumes/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/star_index_${name}_ref/ --genomeFastaFiles $ref_fasta_alt --sjdbGTFfile $ref_gtf --genomeSAindexNbases 13

STAR --runMode genomeGenerate --genomeDir /Volumes/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/star_index_${name}_alt/ --genomeFastaFiles $alt_fasta_alt --sjdbGTFfile $alt_gtf --genomeSAindexNbases 12


#ALIGN THE PAIRED END RNA SEQ READS WITH BOTH THE REFERENCE AND ALTERNATE REFERENCE GENOME
#NOT ALLOWING ANY MULTIMAPPING POSITIONS - SECONDARY ALIGNMENTS WILL STILL BE PRESENT
STAR --genomeDir star_index_${name}_ref/ --readFilesIn $fastq1_base_qual_check $fastq2_base_qual_check --outFilterType BySJout --outFilterMultimapNmax 1 --outFileNamePrefix ${fastq_ID}_${ref_GCA}_filtered_reads_ --outSAMtype BAM SortedByCoordinate && samtools index $ref_aligned_RNA_reads

STAR --genomeDir star_index_${name}_alt/ --readFilesIn $fastq1_base_qual_check $fastq2_base_qual_check --outFilterType BySJout --outFilterMultimapNmax 1 --outFileNamePrefix ${fastq_ID}_${alt_GCA}_filtered_reads_ --outSAMtype BAM SortedByCoordinate && samtools index $alt_aligned_RNA_reads

#FILTERING TO ONLY KEEP TRUE PAIRS
#INPUT
ref_aligned_RNA_reads=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/${fastq_ID}_${ref_GCA}_filtered_reads_Aligned.sortedByCoord.out.bam
alt_aligned_RNA_reads=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/${fastq_ID}_${alt_GCA}_filtered_reads_Aligned.sortedByCoord.out.bam
#OUTPUT
ref_paired_aligned_RNA_reads=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/${fastq_ID}_${ref_GCA}_paired_filtered_reads_Aligned.sortedByCoord.out.bam
alt_paired_aligned_RNA_reads=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/${fastq_ID}_${alt_GCA}_paired_filtered_reads_Aligned.sortedByCoord.out.bam

samtools view -f 2 -b $ref_aligned_RNA_reads > $ref_paired_aligned_RNA_reads && samtools index $ref_paired_aligned_RNA_reads 
samtools view -f 2 -b $alt_aligned_RNA_reads > $alt_paired_aligned_RNA_reads && samtools index $alt_paired_aligned_RNA_reads

#TO LABEL OUR PCR DUPLICATES WE NEED TO USE THIS SCRIPT TO ADD READ GROUP INFORMATION TO THE REF AND ALT ALIGNED RNA SEQ BAM FILES
#INPUT
ref_paired_aligned_RNA_reads=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/${fastq_ID}_${ref_GCA}_paired_filtered_reads_Aligned.sortedByCoord.out.bam
alt_paired_aligned_RNA_reads=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/${fastq_ID}_${alt_GCA}_paired_filtered_reads_Aligned.sortedByCoord.out.bam
#OUTPUT
ref_paired_aligned_RNA_RG_reads=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/${fastq_ID}_${ref_GCA}_paired_filtered_RG_reads_Aligned.sortedByCoord.out.bam
alt_paired_aligned_RNA_RG_reads=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/${fastq_ID}_${alt_GCA}_paired_filtered_RG_reads_Aligned.sortedByCoord.out.bam

picard AddOrReplaceReadGroups I=$ref_paired_aligned_RNA_reads O=$ref_paired_aligned_RNA_RG_reads RGID=sample1 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=sample1 && samtools index $ref_paired_aligned_RNA_RG_reads
picard AddOrReplaceReadGroups I=$alt_aligned_RNA_reads O=$alt_paired_aligned_RNA_RG_reads RGID=sample1 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=sample1 && samtools index $alt_paired_aligned_RNA_RG_reads

#NOW WE CAN ADD FLAGS FOR READS THAT ARE PCR DUPLICATES
#INPUT
ref_paired_aligned_RNA_RG_reads=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/${fastq_ID}_${ref_GCA}_paired_filtered_RG_reads_Aligned.sortedByCoord.out.bam
alt_paired_aligned_RNA_RG_reads=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/${fastq_ID}_${alt_GCA}_paired_filtered_RG_reads_Aligned.sortedByCoord.out.bam
#OUTPUT
ref_aligned_RNA_reads_PCRDUPS=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/${fastq_ID}_${ref_GCA}_paired_filtered_RG_reads_Aligned_duplicates_annotated.sortedByCoord.out.bam
alt_aligned_RNA_reads_PCRDUPS=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/${fastq_ID}_${alt_GCA}_paired_filtered_RG_reads_Aligned_duplicates_annotated.sortedByCoord.out.bam
#TOOL METRICS
picard_ref_metrics=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/Picard_report/${fastq_ID}_${ref_GCA}_dup_metrics.txt
picard_alt_metrics=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/Picard_report/${fastq_ID}_${alt_GCA}_dup_metrics.txt

picard MarkDuplicates I=$ref_paired_aligned_RNA_RG_reads O=$ref_aligned_RNA_reads_PCRDUPS M=$picard_ref_metrics REMOVE_DUPLICATES=false MAX_OPTICAL_DUPLICATE_SET_SIZE=-1 && samtools index $ref_aligned_RNA_reads_PCRDUPS
picard MarkDuplicates I=$alt_paired_aligned_RNA_RG_reads O=$alt_aligned_RNA_reads_PCRDUPS M=$picard_ref_metrics REMOVE_DUPLICATES=false MAX_OPTICAL_DUPLICATE_SET_SIZE=-1 && samtools index $alt_aligned_RNA_reads_PCRDUPS

#GET EXONS THAT ARE ACTUALLY ALIGNED BETWEEN THE TWO HAPLOTYPES
#would be best if this worked with the gft? - i need to pick either gtf or gff 
#INPUT
ref_alt_genome_alignment_bamfile
chromosome_converter_file
ref_annotation
alt_annotation
ref_fasta_index

#OUTPUT
aligned_exons

#SUBSETTING THE BAM FILE OF ALIGNED RNA SEQ READS TO KEEP ONLY READS THAT OVERLAP THE ALIGNED EXONS
#FIRSTLY SPLIT THE ALIGNED EXONS SO THAT WE HAVE REFERENCE AND ALTERNATE SPECIFIC FILES
#(BEDTOOLS INTERSECT NEEDS THEM TO BE SEPERATED)
aligned_exons=/home/s1929681/Desktop/filtered_aligned_exons.csv
ref_aligned_exons_bed=/home/s1929681/Desktop/filtered_aligned_reference_exons.bed
alt_aligned_exons_bed=/home/s1929681/Desktop/filtered_aligned_alternate_exons.bed

sed '1d' $aligned_exons | awk -F',' '{print$2"\t"($4 -1)"\t"$5}' | sort -k1,1n -k2,2n > $ref_aligned_exons_bed 
sed '1d' $aligned_exons | awk -F',' '{print$9"\t"($11 -1)"\t"$12}' | sort -k1,1 -k2,2n > $alt_aligned_exons_bed 

#BEDTOOLS INTERSECT
#INPUTS
ref_aligned_RNA_reads_PCRDUPS=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/${fastq_ID}_${ref_GCA}_paired_filtered_RG_reads_Aligned_duplicates_annotated.sortedByCoord.out.bam
alt_aligned_RNA_reads_PCRDUPS=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/${fastq_ID}_${alt_GCA}_paired_filtered_RG_reads_Aligned_duplicates_annotated.sortedByCoord.out.bam
#OUTPUTS
ref_bamfile_with_only_exon_overlapping_reads=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/${fastq_ID}_${ref_GCA}_Aligned_exon_overlap_paired_filtered_RG_reads_Aligned_duplicates_annotated.sortedByCoord.out.bam
alt_bamfile_with_only_exon_overlapping_reads=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/${fastq_ID}_${alt_GCA}_Aligned_exon_overlap_paired_filtered_RG_reads_Aligned_duplicates_annotated.sortedByCoord.out.bam

bedtools intersect -wa -f 0.75 -a $ref_aligned_RNA_reads_PCRDUPS -b $ref_aligned_exons_bed > $ref_bamfile_with_only_exon_overlapping_reads && samtools index $ref_bamfile_with_only_exon_overlapping_reads

bedtools intersect -wa -f 0.75 -a $alt_aligned_RNA_reads_PCRDUPS -b $alt_aligned_exons_bed > $alt_bamfile_with_only_exon_overlapping_reads && samtools index $alt_bamfile_with_only_exon_overlapping_reads


#PHASE THE READS FOR THE REFERENCE AND ALTERNATE TO READS THAT CAN ACTUALLY BE SPLIT BETWEEN THE REFERENCE AND ALTERNATE
#INPUTS
VCF_ref=/media/s1929681/Seagate/Frankie_DTOL_lep_project/outputs/VCF/${name}.vcf.gz
VCF_alt=/media/s1929681/Seagate/Frankie_DTOL_lep_project/outputs/VCF/${name}_alt_ref.vcf.gz
Chroms=/home/s1929681/One_drive_file_copies_25_01_15/RA_work/Indel_Project/data_set_chroms/${name}_chromosomes.txt

#OUPUTS
#!NEED TO MAKE DEDICATED LOCATION!

python3 /home/s1929681/Desktop/Simon_rna_seq_code/phaseBam.py -b $ref_bamfile_with_only_exon_overlapping_reads -v $VCF_ref -o phase_bam_${name}_${ref_GCA}_filtered_exon_only_reads --use_REF_and_ALT --vcf_to_bam_chrom_file $Chroms 
python3 /home/s1929681/Desktop/Simon_rna_seq_code/phaseBam.py -b $alt_bamfile_with_only_exon_overlapping_reads -v $VCF_alt -o phase_bam_${name}_${alt_GCA}_filtered_exon_only_reads --use_REF_and_ALT 

#WE THEN WHAT TO COUNT THE NUMBER OF READS ALIGNED TO EACH GENE IN THE REFERENCE AND ALT HAP 1 FILES
#question for here is do i need to have created a subsetted version of the gtf file that only has the genes/exons we are looking at - could get aligning exons to output that? - or just filter it in r to only keep the aligned exon genes? 
#INPUT
ref_gtf=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/Bombus_pratorum-${ref_GCA}-2022_03-genes.gtf
alt_gtf=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/Bombus_pratorum-${alt_GCA}-2022_03-genes.gtf
REFERENCE_bamfile_hap_1=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/phase_bam_${name}_${ref_GCA}_filtered_exon_only_reads_1.bam
ALTERNATE_bamfile_hap_1=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/phase_bam_${name}_${alt_GCA}_filtered_exon_only_reads_1.bam

#OUTPUT
ref_feature_counts=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/featureCounts/${name}_${fastq_ID}_${ref_GCA}_pcrdups_aligned_exon_annotated_featurecounts.txt
alt_feature_counts=/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/featureCounts/${name}_${fastq_ID}_${alt_GCA}_pcrdups_aligned_exon_annotated_featurecounts.txt

featureCounts -p -a $ref_gtf -Q 20 -s 0 -o $ref_feature_counts $REFERENCE_bamfile_hap_1
featureCounts -p -a $alt_gtf -Q 20 -s 0 -o $alt_feature_counts $ALTERNATE_bamfile_hap_1

#?I STILL NEED TO MAKE A SUBSETTED FILE FOR THE PCR DUPLICATED READS BUT I DON'T KNOW IF I SHOULD JUST DO THAT FOR THE SUBSETTED READS OR ALL OF THE READS?
