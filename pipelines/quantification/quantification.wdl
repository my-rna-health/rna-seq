version development

import "quant_structs.wdl"
import "quant_sample.wdl"

workflow quantification {
    input {
        Map[String, File] salmon_indexes
        Array[String] samples
        String key = "0a1d74f32382b8a154acacc3a024bdce3709"
        Int extract_threads = 4
        String samples_folder
        Int salmon_threads = 2
        Int bootstraps = 128
    }

    #Headers headers = {"gsm": 0, "series": 1, "run": 2, "path": 3, "organism": 4, "model": 5, "layout": 6, "strategy": 7, "title" : 8, "name": 9, "characteristics": 9}
    scatter(gsm in samples) {

        call get_gsm{
            input: gsm = gsm, key = key
        }

        String gse_folder = samples_folder + "/" + get_gsm.runs[0][1]
        String gsm_folder = gse_folder + "/" + get_gsm.runs[0][0]
        Array[String] headers = get_gsm.headers


       #Array[String] headers = get_gsm.headers
       scatter(run in get_gsm.runs) {
            Array[Pair[String, String]] pairs = zip(headers, run)
            Map[String, String] info = as_map(pairs)
            String layout = run[6] #info["layout"] #run[6]
            Boolean is_paired = (layout != "SINGLE")
            String srr = run[2]#info["run"]
            String sra_folder = gsm_folder + "/" + srr #run[2]

            call download {  input: sra = srr }
            call extract {input: sra = download.out, is_paired = is_paired, threads = extract_threads}
            call fastp { input: reads = extract.out, is_paired = is_paired }
            call copy as copy_report {
             input:
                destination = sra_folder + "/report",
                files = [fastp.report_json, fastp.report_html]
            }

            call copy as copy_cleaned_reads {
             input:
                destination = sra_folder + "/reads",
                files = fastp.reads_cleaned
            }

            String organism = run[4]#info["organism"] #run[4]

            call salmon {
                input:
                    index = salmon_indexes[organism],
                    reads = fastp.reads_cleaned,
                    is_paired = is_paired,
                    threads = salmon_threads,
                    bootstraps = bootstraps,
                    run = srr
            }

            call copy as copy_quant{
            input:
               destination = sra_folder,
               files = [salmon.out]
            }

            File quant_folder = copy_quant.out[0]
            File quant_lib = quant_folder + "/" + "lib_format_counts.json"
            File quant = quant_folder + "/" + "quant.sf"

            #QuantifiedRun quantified_run = {"run": srr, "folder": quant_folder, "quant": quant, "lib": quant_lib, "metainfo": pairs}
            Map[String, File] runs = {"run": srr, "folder": quant_folder, "quant": quant, "lib": quant_lib}
            Map[String, String ] metadata = info
        }

       Array[Map[String, File]] runs_files = runs
       Array[Map[String, String]] runs_metadata = metadata

        #QuantifiedGSM quantified_gsm = {"runs": quantified_run, "metadata": get_gsm.gsm_json}

    }

    output {
        Array[Array[Map[String, File]]] runs_files = runs
        Array[Array[Map[String, String]]] runs_metadata = metadata
        #Array[QuantifiedGSM] quantified_gsms = quantified_gsm
    }
}
