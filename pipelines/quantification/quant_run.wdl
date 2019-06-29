version development

import "extract_run.wdl" as extractor

struct QuantifiedRun {
    String run
    File run_folder
    File quant_folder
    File quant
    File lib
    File genes
    File tx2gene
    Map[String, String] metadata
}

workflow quant_run {
    input {
        String run
        String layout
        File salmon_index
        String folder
        File tx2gene
        Map[String, String] metadata = {"run": run, "layout": layout}

        String key = "0a1d74f32382b8a154acacc3a024bdce3709"
        Int extract_threads = 4
        Int salmon_threads = 2
        Int bootstraps = 128
        Boolean copy_cleaned = false
        String prefix = ""
        Boolean aspera_download = true
    }

    call extractor.extract_run as extract_run{
        input:
            layout = layout,
            run =  run,
            folder = folder,
            copy_cleaned = copy_cleaned,
            extract_threads = extract_threads,
            aspera_download = aspera_download
    }

    call salmon {
        input:
            index = salmon_index,
            reads = extract_run.out.cleaned_reads,
            is_paired = extract_run.out.is_paired,
            threads = salmon_threads,
            bootstraps = bootstraps,
            name = prefix + run
    }

    call tximport {
        input:
            tx2gene =  tx2gene,
            samples = salmon.quant,
            name =  prefix + run
    }

    call extractor.copy as copy_quant{
    input:
       destination = extract_run.out.folder,
       files = [salmon.out, tximport.transcripts, tximport.genes, tximport.genes_length,  tximport.genes_counts ]
    }

    File quant_folder = copy_quant.out[0]
    File quant = quant_folder + "/" + "quant.sf"
    File quant_lib = quant_folder + "/" + "lib_format_counts.json"
    File genes = copy_quant.out[2]

    QuantifiedRun quantified = object {
            run: extract_run.out.run,
            run_folder:extract_run.out.folder,
            quant_folder: quant_folder,
            quant: quant,
            lib: quant_lib,
            genes: genes,
            metadata: metadata,
            tx2gene: tx2gene
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
    File index
    Array[File] reads
    Boolean is_paired
    Int threads
    Int bootstraps = 128
    String name
  }

  command {
    salmon --no-version-check quant -i ~{index}  --numBootstraps ~{bootstraps} --threads ~{threads} -l A --seqBias --gcBias --validateMappings --writeUnmappedNames -o quant_~{name} \
    ~{if(is_paired) then "-1 " + reads[0] + " -2 "+ reads[1] else "-r " + reads[0]}
  }
  # --validateMappings --rangeFactorizationBins ~{rangeFactorizationBins}

  runtime {
    docker: "combinelab/salmon:0.14.1"
    maxRetries: 3
  }

  output {
    File out = "quant_" + name
    File lib = out + "/" + "lib_format_counts.json"
    File quant = out + "/" + "quant.sf"
  }
}


task tximport {
    input {
        File tx2gene
        File samples
        String name
    }

    command {
        /home/rstudio/convert.R --samples ~{samples} --transcripts2genes ~{tx2gene} --name ~{name} --folder expressions
    }

    runtime {
        docker: "quay.io/comp-bio-aging/diff-express:latest"
    }

    output {
        File transcripts = "expressions/transcripts/" + name + "_transcripts_abundance.tsv"

        File genes_length = "expressions/genes/" + name + "_genes_length.tsv"
        File genes_counts = "expressions/genes/" + name + "_genes_counts.tsv"
        File genes = "expressions/genes/" + name + "_genes_abundance.tsv"
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
