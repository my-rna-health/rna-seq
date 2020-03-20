version development

import "quant_run.wdl" as runner
import "extract_run.wdl" as extractor

workflow quant_by_runs{
 input {
        String title = ""
        Array[String] runs
        Map[String, Directory] salmon_indexes
        Map[String, File] transcripts2genes
        String samples_folder
        String key = "0a1d74f32382b8a154acacc3a024bdce3709"
        Int extract_threads = 4
        Int salmon_threads = 4
        Float salmon_max_memory = 20
        Int bootstraps = 96
        Boolean copy_cleaned = false
        Boolean aspera_download = true
    }

    scatter(run in runs) {

        call get_meta{
            input:
                sra = run,
                key = key
        }
        Array[File] metas = get_meta.info

        scatter(json in metas) {

            Map[String, String] info = read_json(json)

            String layout = info["LibraryLayout"]
            Boolean is_paired = (layout != "SINGLE")
            String bioproject = info["BioProject"]
            String experiment = info["Experiment"]
            String organism = info["ScientificName"]
            Directory salmon_index = salmon_indexes[organism]
            File tx2gene = transcripts2genes[organism]

            String sra_folder = samples_folder + "/" + bioproject + "/" + experiment + "/" + run

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
                    metadata = info,
                    gene_map = tx2gene,
                    prefix = bioproject + "_" + experiment + "_",
                    aspera_download = aspera_download
            }
        }

    }

}

task get_meta {
    input{
        String sra
        String key
    }

    command {
       /opt/docker/bin/geo-fetch sra ~{sra} --key ~{key} --output ~{sra}.flat.json
    }

    runtime {
        docker: "quay.io/comp-bio-aging/geo-fetch:0.1.0"
        maxRetries: 1
    }

    output {
        Array[File] info = glob("*.json")
    }
}