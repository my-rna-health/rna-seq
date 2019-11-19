version development

import "quant_sample.wdl" as by_sample

workflow quantification {
    input {
        Map[String, Directory] salmon_indexes
        Map[String, File] transcripts2genes
        Array[String] experiments
        String samples_folder
        String key = "0a1d74f32382b8a154acacc3a024bdce3709"
        Int extract_threads = 4
        Int salmon_threads = 4
        Int bootstraps = 128
        Boolean copy_cleaned = false
        String title = ""
        Boolean aspera_download = true
    }

    scatter(experiment in experiments) {
        call by_sample.quant_sample as quant_sample{
            input:
                experiment = experiment,
                salmon_indexes = salmon_indexes,
                transcripts2genes = transcripts2genes,
                samples_folder = samples_folder,
                key = key,
                extract_threads = extract_threads,
                salmon_threads = salmon_threads,
                bootstraps = bootstraps,
                copy_cleaned = copy_cleaned,
                aspera_download = aspera_download
        }
    }

    output {
        Array[QuantifiedSample] results = quant_sample.sample
    }
}
