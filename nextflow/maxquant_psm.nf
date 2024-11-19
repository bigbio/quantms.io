nextflow.enable.dsl=2

params.msms_file 
params.mzml_dir 
params.output_dir
params.chunksize = 1000000
params.output_prefix_file = "psm"
params.file_num = 20
params.partitions = ""

workflow {
    generateResults(params.msms_file, params.output_dir, params.chunksize, params.output_prefix_file)
    extractInfoFromMzml(generateResults.out, params.mzml_dir, params.output_dir, params.file_num, params.partitions)

}
process generateResults {
    input:
    path msmsFile
    path outputDir
    val chunksize
    val output_prefix_file

    output:
    path "**/*.psm.parquet", emit: 'psm'

    script:
    """
    quantmsioc convert-maxquant-psm --msms_file ${msmsFile} --output_folder ${outputDir} --chunksize ${chunksize} --output_prefix_file ${output_prefix_file}
    """
}

process extractInfoFromMzml {
    input:
    path resultsFile
    path mzmlDir
    path outputDir
    val file_num
    val partitions

    script:
    if (partitions != ''){
        """
        quantmsioc map-spectrum-message-to-parquet --parquet_path ${resultsFile} --mzml_directory ${mzmlDir} --output_folder res/${outputDir} --file_num ${file_num} --partitions ${partitions}
        """
    }else {
        """
        quantmsioc map-spectrum-message-to-parquet --parquet_path ${resultsFile} --mzml_directory ${mzmlDir} --output_folder res/${outputDir} --file_num ${file_num}
        """
    }
}
