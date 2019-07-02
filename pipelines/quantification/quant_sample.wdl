version development

import "quant_run.wdl" as runner
import "extract_run.wdl" as extractor

struct QuantifiedSample {
    String gsm
    String gse
    Array[QuantifiedRun] runs
    File metadata
}

workflow quant_sample {
 input {
        String gsm
        Map[String, File] salmon_indexes
        Map[String, File] transcripts2genes
        String samples_folder
        String key = "0a1d74f32382b8a154acacc3a024bdce3709"
        Int extract_threads = 4
        Int salmon_threads = 2
        Int bootstraps = 128
        Boolean copy_cleaned = false
        Boolean experiment_package = false
        Boolean aspera_download = true
    }

    call get_gsm{
            input: gsm = gsm, key = key, experiment_package = experiment_package
        }


    String gse = sub(get_gsm.runs[0][1], ";", "_")
    String gse_folder = samples_folder + "/" + gse
    String gsm_folder = gse_folder + "/" + gsm
    Array[String] headers = get_gsm.headers

    call extractor.copy as copy_metadata{
    input:
       destination = gsm_folder,
       files = [get_gsm.gsm_json, get_gsm.runs_tsv]
    }

   #Array[String] headers = get_gsm.headers
   scatter(run in get_gsm.runs) {

        Array[Pair[String, String]] pairs = zip(headers, run)
        Map[String, String] info = as_map(pairs)
        String organism = info["organism"] #run[4]
        File salmon_index = salmon_indexes[organism]
        File tx2gene = transcripts2genes[organism]

        String layout = info["layout"] #run[6]
        String srr =  info["run"] #run[2]
        Boolean is_paired = (layout != "SINGLE")
        String sra_folder = gsm_folder + "/" + srr #run[2]


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
                run = srr,
                metadata = info,
                tx2gene = tx2gene,
                prefix = gse + "_" + gsm + "_",
                aspera_download = aspera_download
        }

        Array[Pair[String, String]] pairs = zip(headers, run)
        Map[String, String] info = as_map(pairs)
    }

    output {
        QuantifiedSample sample = object {gsm: gsm, gse: gse, runs: quant_run.quantified_run, metadata: copy_metadata.out[0]}
    }
}

task get_gsm {

    input {
       String gsm
       String key
       Boolean experiment_package = false
    }

    String runs_path = gsm +"_runs.tsv"
    String runs_tail_path = gsm +"_runs_tail.tsv"
    String runs_head_path = gsm +"_runs_head.tsv"


    command {
        /opt/docker/bin/geo-fetch ~{if(experiment_package) then "bioproject --title  N/A" + " --characteristics N/A" + gsm + " " else "gsm"} --key ~{key} -e --output ~{gsm}.json --runs ~{runs_path}  ~{gsm}
        head -n 1 ~{runs_path} > ~{runs_head_path}
        tail -n +2 ~{runs_path} > ~{runs_tail_path}
    }

    runtime {
        docker: "quay.io/comp-bio-aging/geo-fetch:0.0.4"
    }

    output {
        File runs_tsv = runs_path
        File gsm_json = gsm + ".json"
        Array[String] headers = read_tsv(runs_head_path)[0]
        Array[Array[String]] runs = read_tsv(runs_tail_path)
    }
}