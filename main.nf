#!/usr/bin/env nextflow

params.output = "./output/"
params.meta_file = '' // Path to meta.txt
params.count_file = '' // Path to count.txt
params.network_file = ''

params.logFC = true
params.logFC_up = 1
params.logFC_down = -1
params.p_adj = true
params.alpha = 0.05

// data
meta_file = file(params.meta_file)
count_file = file(params.count_file)
network_file = file(params.network_file)

// scripts
kpm_analysis_script = file("${projectDir}/KPMAnalysis.R")


process kpm_analysis {
    container 'kadam0/kpmanalysis:0.0.2'
    publishDir params.output, mode: "copy"

    input:
    path script_file
    path meta_file
    path count_file
    path network_file

    output:                                
    path "*"

    script:
    """
    Rscript $script_file --meta_file ${meta_file} --count_file ${count_file} --network_file ${network_file} --out_dir ${params.output} --logFC ${params.logFC} --logFC_up ${params.logFC_up} --logFC_down ${params.logFC_down} --p_adj ${params.p_adj} --alpha ${params.alpha}
    """
}

workflow {
  kpm_analysis(kpm_analysis_script, meta_file, count_file, network_file)
}
