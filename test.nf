#!/usr/bin/env nextflow

params.reads = "/data/gpfs/assoc/pgl/data/Brassica/Ian_reorder/hic/*_R{1,2}.fastq.gz"
params.genome = "/data/gpfs/assoc/pgl/data/Brassica/Ian_reorder/test/genome.SPLIT.fasta"
params.outdir = "results"

process refIndex {
    conda 'bioconda::chromap=0.2.4'
    input:
        file genome

    output:
        file index

    script:
        """
        ./chromap -i -r $genome -o $index
        """
}

process mapReads {
    tag "$pair_id"

    input:
        file index,
        tuple val(pair_id), path(reads)

    output:
        path pair_id

    script:
        """
        chromap --preset hic -x $index -r ${genome} -1 ${reads[0]} -2 ${reads[1]} --SAM -o $pair_id
        """
    dependsOn refIndex
}

workflow {
    Channel.fromFilePairs(params.reads, checkIfExists: true)
    mapReads(index=Refindex.index, reads=reads, pair_id=reads.basename)
}