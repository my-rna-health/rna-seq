version development
import "https://raw.githubusercontent.com/my-rna-health/rna-seq/master/pipelines/quantification_ena/quant_sample_ena.wdl" as quant


workflow quantification_ena {
    input {
        #Map[String, Directory] salmon_indexes
        #Map[String, File] transcripts2genes
        Directory salmon_index
        File gene_map
        Array[String] samples
        String? email
        String destination
        Int salmon_threads = 4
        Float max_memory = 19
        Int bootstraps = 96
        String prefix = ""
    }

    scatter (sample in samples){
        call quant.quant_sample_ena as quantify {
            input: salmon_index = salmon_index,
                gene_map = gene_map,
                sample = sample,
                email = email,
                destination = destination,
                salmon_threads = salmon_threads,
                max_memory = max_memory,
                bootstraps = bootstraps,
                prefix = prefix

        }
    }

    output {
        Array[Array[MappedRun]] out = quantify.out
    }
}