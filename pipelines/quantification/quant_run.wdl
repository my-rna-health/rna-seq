version development
import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files
import "https://raw.githubusercontent.com/my-rna-health/rna-seq/master/pipelines/quantification/extract_run.wdl" as extractor
#import "extract_run.wdl" as extractor

struct QuantifiedRun {
    String run
    File run_folder
    File quant_folder
    File quant
    File quant_genes
    File gene_map
    File lib
    Map[String, String] metadata
}

workflow quant_run {
    input {
        String run
        String layout
        Directory salmon_index
        String folder
        File gene_map
        Map[String, String] metadata = {"run": run, "layout": layout}

        String key = "0a1d74f32382b8a154acacc3a024bdce3709"
        Int extract_threads = 4
        Int salmon_threads = 4
        Float max_memory = 19
        Int bootstraps = 96
        Boolean copy_cleaned = false
        String prefix = ""
        Boolean aspera_download = true
        Boolean original_names = false
    }

    call extractor.extract_run as extract_run{
        input:
            layout = layout,
            run =  run,
            folder = folder,
            copy_cleaned = copy_cleaned,
            extract_threads = extract_threads,
            aspera_download = aspera_download,
            original_name = original_names
    }

    call salmon {
        input:
            index = salmon_index,
            reads = extract_run.out.cleaned_reads,
            is_paired = extract_run.out.is_paired,
            threads = salmon_threads,
            bootstraps = bootstraps,
            name = prefix + run,
            gene_map = gene_map,
            max_memory = max_memory
    }


    call files.copy as copy_quant{
    input:
       destination = extract_run.out.folder,
       files = [salmon.out]
    }

    File quant_folder = copy_quant.out[0]
    File quant = quant_folder + "/" + "quant.sf"
    File quantified_genes = quant_folder + "/" + "quant.genes.sf"

    File quant_lib = quant_folder + "/" + "lib_format_counts.json"

    QuantifiedRun quantified = object {
            run: extract_run.out.run,
            run_folder: extract_run.out.folder,
            quant_folder: quant_folder,
            quant: quant,
            quant_genes: quantified_genes,
            gene_map: gene_map,
            lib: quant_lib,
            metadata: metadata,
            }

    call write_quant{
        input: quantified_run = quantified, destination = extract_run.out.folder
    }


    output {
        QuantifiedRun quantified_run = quantified
        File quantified_run_json = write_quant.out
    }

}

task salmon {
  input {
    Directory index
    Array[File] reads
    Boolean is_paired
    Int threads
    Int bootstraps = 64
    String name
    File gene_map
    Float max_memory
  }


  command {
    salmon --no-version-check quant -i ~{index} --geneMap ~{gene_map} --numBootstraps ~{bootstraps} --threads ~{threads} -l A --seqBias --gcBias --validateMappings --writeUnmappedNames -o quant_~{name} \
    ~{if(is_paired) then "-1 " + reads[0] + " -2 "+ reads[1] else "-r " + reads[0]}
  }

  runtime {
    docker: "quay.io/biocontainers/salmon@sha256:e56485bfa26913aebaa6351b2ddb1308d0dc0352bf15e7f5431bc58ba5465809" #1.9.0--h7e5ed60_1
    docker_memory: "~{max_memory}G"
    docker_swap: "~{max_memory * 2}G"
    docker_cpu: "~{threads}"
    maxRetries: 3
  }

  output {
    File out = "quant_" + name
    File lib = out + "/" + "lib_format_counts.json"
    File quant = out + "/" + "quant.sf"
    File quant_genes =  out + "/" + "quant.genes.sf"
  }
}


task write_quant {
    input {
        QuantifiedRun quantified_run
        File destination
    }

    String where = sub(destination, ";", "_")

    command {
        cp -L -R -u ~{write_json(quantified_run)} ~{where}/~{quantified_run.run}.json
    }

    output {
        File out = where + "/" + quantified_run.run + ".json"
    }
}