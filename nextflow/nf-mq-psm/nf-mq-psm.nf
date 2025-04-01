// main.nf
nextflow.enable.dsl = 2

// Define default parameters that can be overridden in config
params.raw_dir = params.raw_dir ?: './raw'
params.mzml_dir = params.mzml_dir ?: './mzml'
params.output_dir = params.output_dir ?: './results'
params.msms_file = params.msms_file ?: 'msms.txt'

// Default values, fallback if not specified in config
params.chunksize = params.chunksize ?: 1000000
params.output_prefix_file = params.output_prefix_file ?: 'psm'
params.file_num = params.file_num ?: 20
params.partitions = params.partitions ?: ''

// Log parameters
log.info """\
    MASS SPECTROMETRY WORKFLOW
    =========================
    Raw files directory    : ${params.raw_dir}
    mzML output directory  : ${params.mzml_dir}
    Results directory      : ${params.output_dir}
    MSMS file              : ${params.msms_file}
    Chunk size             : ${params.chunksize}
    Output prefix          : ${params.output_prefix_file}
    File number            : ${params.file_num}
    Partitions             : ${params.partitions ?: 'Not specified'}
    """
    .stripIndent()

process MSCONVERT {
    tag "${rawFile.simpleName}"
    label 'process_medium'
    publishDir "${params.mzml_dir}", mode: 'copy', pattern: '*.mzML', overwrite: true

    container {
        if (workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container) {
            'https://containers.biocontainers.pro/s3/SingImgsRepo/thermorawfileparser/1.3.3--h1341992_0/thermorawfileparser:1.3.3--h1341992_0'
        } else {
            'quay.io/biocontainers/thermorawfileparser:1.3.3--h1341992_0'
        }
    }

    input:
    path rawFile

    output:
    path '*.mzML', emit: mzml

    script:
    """
    ThermoRawFileParser.sh -i=${rawFile} -f=2 -o=./
    """
}

process GENERATE_RESULTS {
    tag "Processing ${msmsFile.simpleName}"
    label 'process_high'
    publishDir "${params.output_dir}/psm", mode: 'copy', pattern: "**/*.psm.parquet"

    input:
    path msmsFile
    val outputDir

    output:
    path "**/*.psm.parquet", emit: psm

    script:
    """
    quantmsioc convert-maxquant-psm \\
        --msms_file ${msmsFile} \\
        --output_folder ${outputDir} \\
        --chunksize ${params.chunksize} \\
        --output_prefix_file ${params.output_prefix_file}
    """
}

process EXTRACT_INFO_FROM_MZML {
    tag "Extracting from ${resultsFile.simpleName}"
    label 'process_high_memory'
    publishDir "${params.output_dir}/extracted", mode: 'copy'

    input:
    path resultsFile
    path mzmlFiles
    val outputDir

    output:
    path "res/${outputDir}/**", emit: extracted_data, optional: true

    script:
    def partitions_param = params.partitions ? "--partitions ${params.partitions}" : ""
    """
    mkdir -p res/${outputDir}
    quantmsioc map-spectrum-message-to-parquet \\
        --parquet_path ${resultsFile} \\
        --mzml_directory ./ \\
        --output_folder res/${outputDir} \\
        --file_num ${params.file_num} \\
        ${partitions_param}
    """
}

workflow {
    // Input channel for raw files
    Channel
        .fromPath("${params.raw_dir}/*.raw")
        .ifEmpty { error "No RAW files found in ${params.raw_dir}" }
        .set { raw_files_ch }

    // Process RAW files to mzML
    mzml_files_ch = MSCONVERT(raw_files_ch)

    // Generate results from msms file
    psm_results_ch = GENERATE_RESULTS(params.msms_file, params.output_dir)

    // Extract information from mzML files using the PSM results
    extracted_data_ch = EXTRACT_INFO_FROM_MZML(
        psm_results_ch,
        mzml_files_ch.collect(),
        params.output_dir
    )
}

// Completion message
workflow.onComplete {
    log.info "Pipeline completed at: ${workflow.complete}"
    log.info "Execution status: ${workflow.success ? 'SUCCESS' : 'FAILED'}"
    log.info "Execution duration: ${workflow.duration}"
}