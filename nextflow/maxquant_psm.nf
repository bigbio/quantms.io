nextflow.enable.dsl=2


process msconvert {
    publishDir "${params.mzml_dir}", mode:'copy', overwrite: true

    if (workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container) {
        container 'https://containers.biocontainers.pro/s3/SingImgsRepo/thermorawfileparser/1.3.3--h1341992_0/thermorawfileparser:1.3.3--h1341992_0'
    }
    else {
        container 'quay.io/biocontainers/thermorawfileparser:1.3.3--h1341992_0'
    }

    input:
    path rawFile

    output:
    path '*.mzML', emit: mzmlFiles

    script:
    """
    mkdir mzml
    ThermoRawFileParser.sh -i=${rawFile} -f=2 -o=./
    """
}

process generateResults {
    input:
    path msmsFile
    path outputDir
    output:
    path "**/*.psm.parquet", emit: 'psm'

    script:
    """
    quantmsioc convert-maxquant-psm --msms_file ${msmsFile} --output_folder ${outputDir} --chunksize ${params.chunksize} --output_prefix_file ${params.output_prefix_file}
    """
}

process extractInfoFromMzml {
    input:
    path resultsFile
    path mzmlDir
    path outputDir

    script:
    if (params.partitions){
        """
        quantmsioc map-spectrum-message-to-parquet --parquet_path ${resultsFile} --mzml_directory ./ --output_folder res/${outputDir} --file_num ${params.file_num} --partitions ${params.partitions}
        """
    }else {
        """
        quantmsioc map-spectrum-message-to-parquet --parquet_path ${resultsFile} --mzml_directory ./ --output_folder res/${outputDir} --file_num ${params.file_num}
        """
    }
}

workflow {
    Channel
    .fromPath("${params.raw_dir}/*.raw")
    .set{ rawFiles }

    msconvert(rawFiles)
    generateResults(params.msms_file, params.output_dir)
    extractInfoFromMzml(generateResults.out, msconvert.out, params.output_dir)

}