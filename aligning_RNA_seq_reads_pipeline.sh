#!/usr/bin/env bash
set -euo pipefail

########################################
# CONFIG
########################################

THREADS=5

SCRATCH_BASE="/scratch/fswift/RNA_seq_reads_api"
OUTPUT_BASE="/mnt/loki/martin/frankie/RNA_seq_glm/1_analysis_files/11_RNA_seq_read_bams"

STAR_REF_BASE="/mnt/loki/martin/frankie/RNA_seq_glm/1_analysis_files/2_STAR_reference_genomes"

########################################
# INPUTS
########################################

name="$1"
RNA_run_accession="$2"
ref_GCA="$3"
alt_GCA="$4"

WORKDIR="${SCRATCH_BASE}/${name}"
FASTQC_DIR="${WORKDIR}/fastQC_reports"
STAR_DIR="${WORKDIR}/STAR_alignment_files"
PICARD_DIR="${WORKDIR}/Picard_report"

mkdir -p "$WORKDIR" "$FASTQC_DIR" "$STAR_DIR" "$PICARD_DIR"

########################################
# LOGGING
########################################

log() {
    echo "[$(date '+%F %T')] $1"
}

########################################
# FUNCTIONS
########################################

download_reads() {
    log "Downloading reads..."

    curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${RNA_run_accession}&result=read_run&fields=fastq_ftp&format=tsv" \
    | tail -n +2 | tr '\t' '\n' | tr ';' '\n' | grep "^ftp" \
    | wget -i - -P "$WORKDIR"
}

rename_reads() {
    log "Renaming reads..."
    cd "$WORKDIR"

    for f in ${RNA_run_accession}*; do
        mv "$f" "${name}_${f}"
    done
}

unzip_reads() {
    log "Unzipping reads..."
    gunzip "${WORKDIR}/${name}_${RNA_run_accession}_1.fastq.gz"
    gunzip "${WORKDIR}/${name}_${RNA_run_accession}_2.fastq.gz"
}

run_fastqc_raw() {
    log "Running FastQC (raw)..."
    fastqc --nogroup \
        "${WORKDIR}/${name}_${RNA_run_accession}_1.fastq" \
        "${WORKDIR}/${name}_${RNA_run_accession}_2.fastq"

    mv "${WORKDIR}"/*fastqc* "$FASTQC_DIR"/
}

run_fastp() {
    log "Running fastp..."

    fastp \
      -i "${WORKDIR}/${name}_${RNA_run_accession}_1.fastq" \
      -I "${WORKDIR}/${name}_${RNA_run_accession}_2.fastq" \
      -o "${WORKDIR}/${RNA_run_accession}_1_base_quaility_checked.fastq" \
      -O "${WORKDIR}/${RNA_run_accession}_2_base_quaility_checked.fastq" \
      -q 20 -l 30 \
      -h "${FASTQC_DIR}/${name}_${RNA_run_accession}_fastp_report.html" \
      -j "${FASTQC_DIR}/${name}_${RNA_run_accession}_fastp_report.json"

    rm -f "${WORKDIR}/${name}_${RNA_run_accession}_1.fastq"
    rm -f "${WORKDIR}/${name}_${RNA_run_accession}_2.fastq"
}

run_fastqc_trimmed() {
    log "Running FastQC (trimmed)..."

    fastqc --nogroup \
        "${WORKDIR}/${RNA_run_accession}_1_base_quaility_checked.fastq" \
        "${WORKDIR}/${RNA_run_accession}_2_base_quaility_checked.fastq"

    mv "${WORKDIR}"/*fastqc* "$FASTQC_DIR"/
}

run_star() {
    local genome_type=$1
    local gca=$2

    log "Running STAR (${genome_type})..."

    STAR \
      --runThreadN "$THREADS" \
      --genomeLoad NoSharedMemory \
      --genomeDir "${STAR_REF_BASE}/star_index_${name}_${genome_type}" \
      --readFilesIn \
        "${WORKDIR}/${RNA_run_accession}_1_base_quaility_checked.fastq" \
        "${WORKDIR}/${RNA_run_accession}_2_base_quaility_checked.fastq" \
      --outFilterType BySJout \
      --outFilterMultimapNmax 1 \
      --outFileNamePrefix "${WORKDIR}/${RNA_run_accession}_${gca}_filtered_reads_" \
      --outSAMmultNmax 1 \
      --outSAMtype BAM SortedByCoordinate

    samtools index "${WORKDIR}/${RNA_run_accession}_${gca}_filtered_reads_Aligned.sortedByCoord.out.bam"
}

filter_proper_pairs() {
    local gca=$1

    log "Filtering proper pairs (${gca})..."

    samtools view -f 2 -b \
      "${WORKDIR}/${RNA_run_accession}_${gca}_filtered_reads_Aligned.sortedByCoord.out.bam" \
      > "${WORKDIR}/${RNA_run_accession}_${gca}_paired_filtered_reads_Aligned.sortedByCoord.out.bam"

    samtools index "${WORKDIR}/${RNA_run_accession}_${gca}_paired_filtered_reads_Aligned.sortedByCoord.out.bam"
}

add_read_groups() {
    local gca=$1

    log "Adding read groups (${gca})..."

    picard AddOrReplaceReadGroups \
      I="${WORKDIR}/${RNA_run_accession}_${gca}_paired_filtered_reads_Aligned.sortedByCoord.out.bam" \
      O="${WORKDIR}/${RNA_run_accession}_${gca}_paired_filtered_RG_reads_Aligned.sortedByCoord.out.bam" \
      RGID=sample1 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=sample1

    samtools index "${WORKDIR}/${RNA_run_accession}_${gca}_paired_filtered_RG_reads_Aligned.sortedByCoord.out.bam"
}

mark_duplicates() {
    local gca=$1

    log "Marking duplicates (${gca})..."

    picard MarkDuplicates \
      I="${WORKDIR}/${RNA_run_accession}_${gca}_paired_filtered_RG_reads_Aligned.sortedByCoord.out.bam" \
      O="${WORKDIR}/${RNA_run_accession}_${gca}_paired_filtered_RG_reads_Aligned_duplicates_annotated.sortedByCoord.out.bam" \
      M="${PICARD_DIR}/${RNA_run_accession}_${gca}_dup_metrics.txt" \
      REMOVE_DUPLICATES=false MAX_OPTICAL_DUPLICATE_SET_SIZE=-1

    samtools index "${WORKDIR}/${RNA_run_accession}_${gca}_paired_filtered_RG_reads_Aligned_duplicates_annotated.sortedByCoord.out.bam"
}

cleanup_intermediate() {
    log "Cleaning up intermediate files..."

    rm -f "$WORKDIR"/*_base_quaility_checked.fastq
    rm -f "$WORKDIR"/*filtered_reads_Aligned.sortedByCoord.out.bam*
    rm -f "$WORKDIR"/*_paired_filtered_reads_Aligned.sortedByCoord.out.bam*
    rm -f "$WORKDIR"/*_paired_filtered_RG_reads_Aligned.sortedByCoord.out.bam*
}

finalise_outputs() {
    log "Finalising outputs..."

    mv "$WORKDIR"/*out "$STAR_DIR"/ 2>/dev/null || true
    mv "$WORKDIR"/*out.tab "$STAR_DIR"/ 2>/dev/null || true

    cd "$WORKDIR"

    log "Renaming BAM files to remove RNA accession..."

    # Rename BAM + BAI files
    for f in ${RNA_run_accession}_*.bam*; do
        new_name="${name}_${f#${RNA_run_accession}_}"
        mv "$f" "$new_name"
    done

    rsync -av "$WORKDIR"/*.bam* "${OUTPUT_BASE}/0_bams/"
    rm -f "$WORKDIR"/*.bam* 
    rsync -av "$WORKDIR" "${OUTPUT_BASE}/1_alignment_extra_files/"
    rm -rf "${STAR_REF_BASE}/star_index_${name}_*"
    rm -rf "$WORKDIR"
}

########################################
# PIPELINE
########################################

main() {
    download_reads
    rename_reads
    unzip_reads

    run_fastqc_raw
    run_fastp
    run_fastqc_trimmed

    run_star "ref" "$ref_GCA"
    run_star "alt" "$alt_GCA"

    filter_proper_pairs "$ref_GCA"
    filter_proper_pairs "$alt_GCA"

    add_read_groups "$ref_GCA"
    add_read_groups "$alt_GCA"

    mark_duplicates "$ref_GCA"
    mark_duplicates "$alt_GCA"

    cleanup_intermediate
    finalise_outputs

    log "Pipeline completed successfully 🎉"
}

main