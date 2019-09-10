version development

import "star_align_run.wdl" as runner
import "extract_run.wdl" as extractor

workflow star_align_by_runs{
 input {
        Array[String] runs
        Map[String, Directory] indexes
        Map[String, File] transcripts2genes
        String samples_folder

        String key = "0a1d74f32382b8a154acacc3a024bdce3709"
        Int extract_threads = 4
        Int salmon_threads = 4
        Int bootstraps = 128
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
            String organism = info["ScientificName"]
            Directory index = indexes[organism]
            File tx2gene = transcripts2genes[organism]

            String sra_folder = samples_folder + "/" + "bioprojects" + "/" + bioproject + "/" + run

            call runner.star_align_run as star_align_run{
                input:
                    bootstraps = bootstraps,
                    index_dir = index,
                    key = key,
                    layout = layout,
                    folder = sra_folder,
                    copy_cleaned = copy_cleaned,
                    salmon_threads = salmon_threads,
                    extract_threads = extract_threads,
                    run = run,
                    metadata = info,
                    tx2gene = tx2gene,
                    prefix = bioproject + "_" + run + "_",
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
        docker: "quay.io/comp-bio-aging/geo-fetch:0.0.5"
    }

    output {
        Array[File] info = glob("*.json")
    }
}