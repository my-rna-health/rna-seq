version development
import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files
import "https://raw.githubusercontent.com/my-rna-health/rna-seq/master/pipelines/quantification_ena/download_sample.wdl" as download

struct MappedRun {
    RunInfo run
    File quant_genes
    File quant
    File lib
}

workflow quant_sample_ena{
    input {
        String sample
        String? email
        String destination
        Directory salmon_index
        File gene_map
        Int salmon_threads = 4
        Float max_memory = 19
        Int bootstraps = 96
        String prefix = ""
    }

    call download.download_sample as get_sample{
        input: sample = sample, email = email, destination = destination
    }

    scatter(run in get_sample.runs) {
        call fastp {
            input:
                reads = run.reads,
                is_paired = true
        }

        call salmon {
            input:
                index = salmon_index,
                reads = fastp.reads_cleaned,
                is_paired = true,
                threads = salmon_threads,
                bootstraps = bootstraps,
                name = prefix + run.run_accession,
                gene_map = gene_map,
                max_memory = max_memory
        }

        call files.copy as copy_quant{
            input:
                destination = run.run_folder,
                files = [salmon.out]
        }
        #MappedRun mapped_run = object {
        #    run: run, quant_genes: salmon.quant_genes, quant: salmon.quant, lib: salmon.lib
        #}
    }

    output {
        #Array[MappedRun] mapped_runs = mapped_run
    }
}


task fastp {
    input {
        Array[File] reads
        Boolean is_paired
    }

    command {
        fastp --cut_front --cut_tail --cut_right --trim_poly_g --trim_poly_x --overrepresentation_analysis \
        -i ~{reads[0]} -o ~{basename(reads[0], ".fastq.gz")}_cleaned.fastq.gz \
        ~{if( is_paired ) then "--detect_adapter_for_pe " + "--correction -I "+reads[1]+" -O " + basename(reads[1], ".fastq.gz") +"_cleaned.fastq.gz" else ""}
    }

    runtime {
        docker: "quay.io/biocontainers/fastp@sha256:2489fe56260bde05bdf72a8ead4892033b9a05dc4525affb909405bea7839d1b" #0.23.2--h5f740d0_3
    }

    output {
        File report_json = "fastp.json"
        File report_html = "fastp.html"
        Array[File] reads_cleaned = if( is_paired )
                                    then [basename(reads[0], ".fastq.gz") + "_cleaned.fastq.gz", basename(reads[1], ".fastq.gz") + "_cleaned.fastq.gz"]
                                    else [basename(reads[0], ".fastq.gz") + "_cleaned.fastq.gz"]
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