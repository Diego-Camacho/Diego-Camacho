nextflow.enable.dsl=2

params.input_dir = "."
params.debug_dir = "."
params.script_dir = "."
params.folder = "."

process ALIGN {
    input:
    tuple val(name), path(reads)

    output:
    tuple val(name), path("${name}_aligned.bam"), emit: aligned

    script:
    """
    ${params.script_dir}/2_2_2_p_alignment.sh ${params.folder} ${name}
    """
}

process DELETE_CHRM {
    input:
    tuple val(name), path(bam)

    output:
    tuple val(name), path("${name}_nochrM.bam"), emit: nochrM

    script:
    """
    ${params.script_dir}/2_2_3_p_delchrM.sh ${params.folder} ${name}
    """
}

process MARK_DUPLICATES {
    input:
    tuple val(name), path(bam)

    output:
    tuple val(name), path("${name}_markdup.bam"), emit: markdup

    script:
    """
    ${params.script_dir}/2_2_4_p_picard_markdup.sh ${params.folder} ${name}
    """
}

process INDEX_BAM {
    input:
    tuple val(name), path(bam)

    output:
    tuple val(name), path("${name}_indexed.bam"), emit: indexed

    script:
    """
    ${params.script_dir}/2_2_5_p_index_bam.sh ${params.folder} ${name}
    """
}

process REMOVE_NON_UNIQUE {
    input:
    tuple val(name), path(bam)

    output:
    tuple val(name), path("${name}_unique.bam"), emit: unique

    script:
    """
    ${params.script_dir}/2_2_6_p_remove_NonUniqueAligned.sh ${params.folder} ${name}
    """
}

process SAM_TO_BED {
    input:
    tuple val(name), path(bam)

    output:
    tuple val(name), path("${name}.bed"), emit: bed

    script:
    """
    ${params.script_dir}/2_2_7_p_samtobed.sh ${params.folder} ${name}
    """
}

process MACS2 {
    input:
    tuple val(name), path(bed)

    output:
    tuple val(name), path("${name}_peaks.narrowPeak"), emit: peaks

    script:
    """
    ${params.script_dir}/2_2_8_p_macs2.sh ${params.folder} ${name}
    """
}

workflow {
    Channel
        .fromPath("${params.input_dir}/B*_R1.fastq.gz")
        .map { file -> 
            def name = file.name.toString().tokenize('_')[0]
            tuple(name, file)
        }
        .set { input_files }

    ALIGN(input_files)
    DELETE_CHRM(ALIGN.out.aligned)
    MARK_DUPLICATES(DELETE_CHRM.out.nochrM)
    INDEX_BAM(MARK_DUPLICATES.out.markdup)
    REMOVE_NON_UNIQUE(INDEX_BAM.out.indexed)
    SAM_TO_BED(REMOVE_NON_UNIQUE.out.unique)
    MACS2(SAM_TO_BED.out.bed)
}
