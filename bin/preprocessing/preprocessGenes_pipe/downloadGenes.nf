nextflow.enable.dsl=2

// Processes

// process downloadHumanPC { 

//     publishDir "data/PCgenes", mode: 'symlink'

//     output:
//         path "refseq_pcoding_human.tsv"
//     shell:
//         '''
//         downloadP.R
//         '''
// }

// process matchHGNCSymbol {
//     publishDir "data/PCgenes", mode: 'symlink'

//     input:
//         path refseqPC

//     output:
//         path 'humanPCMaster.tsv'
//     shell:
//         '''
//         getHGNC.R !{refseqPC}
//         '''
// }

process downloadAllGenes { 

    publishDir "data/genes", mode: 'copy'

    output:
        path "humanAllGenes.csv"
        path "humanPCGenes.csv"
    shell:
        '''
        getAllGenes.R
        '''
}

// workflow downloadPCGenes_pipe {
//     take:

//     main:
//         downloadHumanPC()
//         matchHGNCSymbol(downloadHumanPC.out)

//     emit:
//         masterPC = matchHGNCSymbol.out
// }


workflow downloadGenes_pipe {
    take:

    main:
        downloadAllGenes()

    // emit:
    //     masterPC = downloadAllGenes.out
}