version development

import "quant_run.wdl" as runner
import "extract_run.wdl" as extractor

workflow quant_by_runs{
 input {
        Array[String] runs
        Map[String, File] salmon_indexes
        Map[String, File] transcripts2genes
        String samples_folder

        String key = "0a1d74f32382b8a154acacc3a024bdce3709"
        Int extract_threads = 4
        Int salmon_threads = 2
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
        Map[String, String] info = get_meta.info

        String layout = info["LibraryLayout"]
        Boolean is_paired = (layout != "SINGLE")
        String bioproject = info["BioProject"]
        String organism = info["ScientificName"]
        File salmon_index = salmon_indexes[organism]
        File tx2gene = transcripts2genes[organism]

        String sra_folder = samples_folder + "/" + "bioprojects" + "/" + bioproject + "/" + run

        call runner.quant_run as quant_run{
            input:
                bootstraps = bootstraps,
                salmon_index = salmon_index,
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

task get_meta {
    input{
        String sra
        String key
    }

    String runs_path = sra +".tsv"
    String runs_tail_path = sra +"_tail.tsv"
    String runs_head_path = sra +"_head.tsv"

    command {
       /opt/docker/bin/geo-fetch sra ~{sra} --key ~{key} --output ~{sra}.tsv
       head -n 1 ~{runs_path} > ~{runs_head_path}
       tail -n +2 ~{runs_path} > ~{runs_tail_path}
    }

    runtime {
        docker: "quay.io/comp-bio-aging/geo-fetch:0.0.4"
    }

    output {
        File run_tsv = runs_path
        Array[String] headers = read_tsv(runs_head_path)[0]
        Array[Array[String]] run = read_tsv(runs_tail_path)
        Array[Pair[String, String]] pairs = zip(headers, run[0])
        Map[String, String] info = as_map(pairs)
    }
}