version development

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files
import "https://raw.githubusercontent.com/my-rna-health/rna-seq/master/pipelines/quantification/extract_run.wdl" as extractor
#import "download_sample.wdl" as extractor
import "https://raw.githubusercontent.com/my-rna-health/rna-seq/master/pipelines/quantification/quant_run.wdl" as runner
#import "quant_run.wdl" as runner

workflow quant_by_runs_simple{
    input {
        String title = ""
        Array[String] runs
        Map[String, Directory] salmon_indexes
        Map[String, File] transcripts2genes
        String species

        String results_folder

        String key = "0a1d74f32382b8a154acacc3a024bdce3709"
        Int extract_threads = 4
        Int salmon_threads = 4
        Float salmon_max_memory = 20
        Int bootstraps = 96
        Boolean copy_cleaned = false
        Boolean aspera_download = true
        String layout = "PAIRED"
    }

    scatter(run in runs) {
            Boolean is_paired = (layout != "SINGLE")
            Directory salmon_index = salmon_indexes[species]
            File tx2gene = transcripts2genes[species]

            String sra_folder = results_folder + "/" + run

            call runner.quant_run as quant_run{
                input:
                    bootstraps = bootstraps,
                    salmon_index = salmon_index,
                    key = key,
                    layout = layout,
                    folder = sra_folder,
                    copy_cleaned = copy_cleaned,
                    salmon_threads = salmon_threads,
                    max_memory = salmon_max_memory,
                    extract_threads = extract_threads,
                    run = run,
                    gene_map = tx2gene,
                    aspera_download = aspera_download
            }

    }

}
