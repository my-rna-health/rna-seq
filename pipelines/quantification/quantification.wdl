version development

import "quant_sample.wdl" as by_sample

workflow quantification {
    input {
        Map[String, Directory] salmon_indexes
        Map[String, File] transcripts2genes
        Array[String] samples
        String samples_folder
        String key = "0a1d74f32382b8a154acacc3a024bdce3709"
        Int extract_threads = 4
        Int salmon_threads = 4
        Int bootstraps = 128
        Boolean copy_cleaned = false
        String title = ""
    }

    scatter(gsm in samples) {
        call by_sample.quant_sample as quant_sample{
            input:
                gsm = gsm,
                salmon_indexes = salmon_indexes,
                transcripts2genes = transcripts2genes,
                samples_folder = samples_folder,
                key = key,
                extract_threads = extract_threads,
                salmon_threads = salmon_threads,
                bootstraps = bootstraps,
                copy_cleaned = copy_cleaned
        }
    }

    output {
        Array[QuantifiedSample] results = quant_sample.sample
    }
}
