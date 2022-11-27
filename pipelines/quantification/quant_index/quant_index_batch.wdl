version development

#import "quant_index.wdl" as index
import "https://raw.githubusercontent.com/my-rna-health/rna-seq/master/pipelines/quantification/quant_index/quant_index.wdl" as index

workflow quant_index_batch{
    input {
        String indexes_folder
        Array[Array[Gentrome]] references
        Int threads_per_index = 3
        String memory_per_index = "24G"
    }

    scatter(ref in references) {
        call index.quant_index{
            input:
                references = ref,
                indexes_folder = indexes_folder,
                threads = threads_per_index,
                max_memory = memory_per_index
        }
    }
}