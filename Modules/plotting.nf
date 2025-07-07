nextflow.enable.dsl=2

process plotPKS {

    label 'process_low'
    scratch true
    publishDir "${params.pks_dir}", mode: 'copy'
    conda "${params.pks_env}"

    input:
    path coverage_file

    output:
    path "*.pks.cirocs.pdf"

    script:
    """
    sample_name=\$(basename ${coverage_file} | cut -d '.' -f1)
    output_pdf="\${sample_name}.pks.cirocs.pdf"

    Rscript "${params.scripts}/plotPKS.R" ${coverage_file} "${params.pks_cytoband}" "\${output_pdf}"
    """
}

process plotECOLI {

    label 'process_low'
    scratch true
    publishDir "${params.pks_dir}", mode: 'copy'
    conda "${params.pks_env}"

    input:
    path coverage_file

    output:
    path "*.ecoli.cirocs.pdf"

    script:
    """
    sample_name=\$(basename ${coverage_file} | cut -d '.' -f1)
    output_pdf="\${sample_name}.ecoli.cirocs.pdf"

    Rscript "${params.scripts}/plotECOLI.R" ${coverage_file} "${params.ecoli_cytoband}" "\${output_pdf}"
    """
}



process masterTable {
    label 'process_low'
    scratch true
    publishDir "${params.pks_dir}", mode: 'copy'
    conda "${params.pks_env}"

    input:
    path coverage_files

    output:
    path "ecoli.pks.gene.counts.txt"

    script:
    def coverage_str = coverage_files.collect { it.getName() }.join(' ')

    """
    python "${params.scripts}/mergeGeneCounts.py" ${coverage_str} ecoli.pks.gene.counts.txt
    """
}


