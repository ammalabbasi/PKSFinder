nextflow.enable.dsl=2

process pksFinder {
    label 'process_medium'
    scratch true
    publishDir "${params.pks_dir}", mode: 'copy'
    conda "${params.pks_env}"  // Set conda environment

    input:
    tuple path(r1_fastq), path(r2_fastq)
 
    output:
    tuple path("*.coverage.txt"),
          path("*.coverage.bedgraph"),
          path("*.counts.txt"),
          path("*.sorted.bam"),
          path("*.sam")


    script:

    def sample_name = r1_fastq.baseName.split('\\.')[0]
    def bedtools_cov = "${sample_name}.coverage.txt"
    def coverage = "${sample_name}.coverage.bedgraph"
    def counts = "${sample_name}.counts.txt"
    def bam = "${sample_name}.sorted.bam"
    def sam = "${sample_name}.sam"

    """
    if [[ -f "${params.pks_dir}/${bedtools_cov}" && -f "${params.pks_dir}/${coverage}" && -f "${params.pks_dir}/${counts}" && -f "${params.pks_dir}/${bam}" && -f "${params.pks_dir}/${sam}" ]]; then
        echo "Skipping PKSFinder process: Found required files"
        for f in "${bedtools_cov}" "${coverage}" "${counts}" "${bam}" "${sam}"; do
            if [[ ! -f "\$f" ]]; then
                ln -s "${params.pks_dir}/\$f" . 2>/dev/null || cp "${params.pks_dir}/\$f" .
            fi
        done
        exit 0
    fi

    # Actual execution
    # Combine R1 and R2 fastq files into a single file
    cat "${r1_fastq}" "${r2_fastq}" > "${sample_name}.trimmed.fastq.gz"

    echo 'Bowtie2 Alignment Classified'
    bowtie2 -x "${params.pks_genome}" -q "${sample_name}.trimmed.fastq.gz" --seed 42 --threads 1 --very-sensitive --no-unal -S ${sam}
    samtools view -bS -q 40 ${sam} | samtools sort -o ${bam} -
    samtools index ${bam}
    "${params.scripts}/featureCounts" -a "${params.pks_genome_annotation}" -o ${counts} ${bam} -t 'gene' -F 'GTF' -g 'Name'
    bamCoverage -b ${bam} -o ${coverage} --normalizeUsing RPKM --outFileFormat 'bedgraph'
    samtools view -b ${bam} | bedtools genomecov -ibam - -d > ${bedtools_cov}
    
    """
}