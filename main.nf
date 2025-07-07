nextflow.enable.dsl=2

// Edit this with your sample.csv path
params.sample = "/tscc/nfs/home/amabbasi/restricted/microbiome_pipeline/sample.csv"

// Output directories
params.unmapped_bam_dir = "${projectDir}/RESULTS/UNMAPPED_BAM"
params.mapped_reads_dir = "${projectDir}/RESULTS/MAPPED_READS"
params.pks_dir = "${projectDir}/RESULTS/PKS_READS"

// Databases and ref files [CHANGE THIS]
params.hg38_db="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/human-GRC-db.mmi"
params.t2t_phix_db="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/human-GCA-phix-db.mmi"
params.pangenome_db="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/pangenome_mmi"
params.adapters="${projectDir}/ref/known_adapters.fna"
params.pks_genome="${projectDir}/indices/GCF_000025745.1/GCF_000025745.1_ASM2574v1_genomic"
params.pks_genome_annotation="${projectDir}/indices/GCF_000025745.1/genomic.gff"
params.pks_cytoband="${projectDir}/indices/GCF_000025745.1/genomic_pks.txt"
params.ecoli_cytoband="${projectDir}/indices/GCF_000025745.1/genomic_ecoli.txt"

// Enviroment paths
params.samtools_env = "./conda_envs/samtools_env.yml"
params.fastp_env = "./conda_envs/fastp_env.yml"
params.minimap2_env = "./conda_envs/minimap2_env.yml"
params.pks_env = "./conda_envs/pks_env.yml"

// Package and script paths
params.scripts ="${projectDir}/scripts"


// Include the external processes
include { extractReads } from './Modules/extract_reads.nf'
include { filterReads } from './Modules/filter_reads.nf'
include { mapReads as mapReads } from './Modules/map_reads.nf'
include { pksFinder } from './Modules/pksFinder.nf'
include { plotPKS;plotECOLI;masterTable } from './Modules/plotting.nf'


// Define the workflow
workflow {

    // ------------------- STEP1: HOST DEPLETION ---------------------- //

    // Read and parse the sample sheet
    sample_sheet = nextflow.Channel.fromPath(params.sample)
        .splitCsv(header: true)
        .map { row ->
            row.subMap('patient', 'bam') // Extract relevant metadata
        }

    // Extract reads from BAM files
    extractReads(sample_sheet).set { UNMAPPED_READS }

    UNMAPPED_READS.multiMap { sampleID, r1, r2 -> 
        path_only: tuple(r1, r2)
        whole: tuple(sampleID, r1, r2)
    }
    .set { UNMAPPED_READS_MULTI }
    
    // Filter poor quality reads using fastp
    filterReads(UNMAPPED_READS_MULTI.whole).set { FILTERED_UNMAPPED_READS }

    FILTERED_UNMAPPED_READS.multiMap { sampleID, r1, r2 ->
        path_only: tuple(r1, r2)
        whole: tuple(sampleID, r1, r2)
    }
    .set { FILTERED_UNMAPPED_READS_MULTI }


    // gather the list of pangenome .mmi files
    def mmiFiles = []
    def dir = new File("${params.pangenome_db}")
    dir.eachFileRecurse (groovy.io.FileType.FILES) { file ->
        if (file.name.endsWith('.mmi')) {
            mmiFiles << file
        }
    }

    // Processing for both READSs
    mapReads(FILTERED_UNMAPPED_READS_MULTI.whole, mmiFiles).set { MAPPED_READS }

    MAPPED_READS.multiMap { sampleID, r1Hg38, r1T2T, r1Pan, r2Hg38, r2T2T, r2Pan ->
        Hg38: tuple(r1Hg38, r2Hg38) 
        T2T: tuple(r1T2T, r2T2T) 
        PAN: tuple(r1Pan, r2Pan) 
    }
    .set { MAPPED_READS_MULTI }


    // ------------------- STEP2: MAP TO PKS ISLAND ---------------------- //

    pksFinder(MAPPED_READS_MULTI.PAN).set { PKS_OUT }

    // ------------------- STEP3: PLOT ECOLI & PKS ISLAND ---------------------- //
    
    PKS_OUT.map {it[1]}
    .set { COVERAGE_BEDGRAPH }

    //plotPKS(COVERAGE_BEDGRAPH)
    plotECOLI(COVERAGE_BEDGRAPH)

    // ------------------- STEP4: GENERATE ECOLI AND PKS MASTER TABLE ---------------------- //

    // Extract count files
    PKS_OUT.map { it[2] }
    .collect()
    .set { COUNT_FILES }
    masterTable(COUNT_FILES)

}



