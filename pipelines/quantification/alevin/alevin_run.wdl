version development

import "extract_run.wdl" as extractor

workflow alevin_run {
    input {
        String run
        File umi
        Directory salmon_index
        String folder
        File gene_map
        String protocol = "chromiumV3"
        String key = "0a1d74f32382b8a154acacc3a024bdce3709"
        Int extract_threads = 8
        Int salmon_threads = 8
        Boolean copy_cleaned = false
        String prefix = ""
        Boolean aspera_download = true
    }

    call extractor.extract_run as extract_run{
        input:
            layout = "SINGLE",
            run =  run,
            folder = folder,
            copy_cleaned = copy_cleaned,
            extract_threads = extract_threads,
            aspera_download = aspera_download
    }

    call alevin {
        input:
            index = salmon_index,
            umi = umi,
            read = extract_run.out.cleaned_reads[0],
            threads = salmon_threads,
            name = prefix + run,
            gene_map = gene_map
    }

    call extractor.copy as copy_quant{
    input:
       destination = extract_run.out.folder,
       files = [alevin.out]
    }

}

task alevin {
  input {
    Directory index
    File umi
    File read
    Boolean is_paired
    Int threads
    String name
    File gene_map
    String protocol
  }


  command {
    salmon alevin -i ~{index} --~{protocol} --tgMap ~{gene_map} --threads ~{threads} -l A -o quant_~{name} \
    ~{"-1 " + umi + " -2 "+ read}
  }

  runtime {
    docker: "quay.io/biocontainers/salmon@sha256:b1b5136321e8d5849e49035cd59e5dda755ba759f4e6fe3ffe1e914444a711af" #1.3.0--hf69c8f4_0  #1.0.0--hf69c8f4_0
    maxRetries: 3
  }

  output {
    File out = "quant_" + name
  }
}