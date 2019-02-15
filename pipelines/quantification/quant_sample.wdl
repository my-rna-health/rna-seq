version development
import "quant_structs.wdl"

workflow quant_sample {
 input {

        String gsm
        Map[String, File] salmon_indexes
        String samples_folder
        String key = "0a1d74f32382b8a154acacc3a024bdce3709"
        Int extract_threads = 4
        Int salmon_threads = 2
        Int bootstraps = 128
    }

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
}

task get_gsm {

    input {
       String gsm
       String key
    }

    String runs_path = gsm +"_runs.tsv"
    String runs_tail_path = gsm +"_runs_tail.tsv"
    String runs_head_path = gsm +"_runs_head.tsv"


    command {
        /opt/docker/bin/geo-fetch gsm --key ~{key} -e --output ~{gsm}.json --runs ~{runs_path}  ~{gsm}
        head -n 1 ~{runs_path} > ~{runs_head_path}
        tail -n +2 ~{runs_path} > ~{runs_tail_path}
    }

    runtime {
        docker: "quay.io/comp-bio-aging/geo-fetch:0.0.2"
    }

    output {
        File runs_tsv = runs_path
        File gsm_json = gsm + ".json"
        Array[String] headers = read_tsv(runs_head_path)[0]
        Array[Array[String]] runs = read_tsv(runs_tail_path)
    }
}

task download {
    input {
        String sra
    }

    command {
        download_sra_aspera.sh ~{sra}
    }

    #https://github.com/antonkulaga/biocontainers/tree/master/downloaders/sra

    runtime {
        docker: "quay.io/antonkulaga/download_sra:latest"
        #maxRetries: 2
    }

    output {
        File out = "results" + "/" + sra + ".sra"
     }
}

task extract {
    input {
        File sra
        Boolean is_paired
        Int threads
    }

    String name = basename(sra, ".sra")
    String folder = "extracted"
    String prefix = folder + "/" + name
    String prefix_sra = prefix + ".sra"

    #see https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump for docs

    command {
        fasterq-dump --outdir ~{folder} --threads ~{threads} --progress --split-files --skip-technical ~{sra}
        ~{if(is_paired) then "mv" + " " + prefix_sra + "_1.fastq" + " " + prefix + "_1.fastq"  else "mv" + " " + prefix_sra + ".fastq" + " " + prefix + ".fastq"}
        ~{if(is_paired) then "mv" + " " + prefix_sra + "_2.fastq" + " " + prefix + "_2.fastq"  else ""}
    }

    runtime {
        docker: "quay.io/biocontainers/sra-tools@sha256:b03fd02fefc3e435cd36eef802cc43decba5d13612142e9bc9610f2727364f4f" #2.9.1_1--h470a237_0
        #maxRetries: 3
    }

    output {
        Array[File] out = if(is_paired) then [prefix + "_1.fastq",  prefix + "_2.fastq"] else [prefix + ".fastq"]
     }
}



task fastp {
    input {
        Array[File] reads
        Boolean is_paired
    }

    command {
        fastp --cut_by_quality5 --cut_by_quality3 --trim_poly_g --overrepresentation_analysis \
            -i ~{reads[0]} -o ~{basename(reads[0], ".fastq.gz")}_cleaned.fastq.gz \
            ~{if( is_paired ) then "--detect_adapter_for_pe " + "--correction -I "+reads[1]+" -O " + basename(reads[1], ".fastq.gz") +"_cleaned.fastq.gz" else ""}
    }

    runtime {
        docker: "quay.io/biocontainers/fastp@sha256:159da35f3a61f6b16650ceef6583c49d73396bc2310c44807a0d929c035d1011" #0.19.5--hd28b015_0
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
    File index
    Array[File] reads
    Boolean is_paired
    Int threads
    Int bootstraps = 128
    String run
  }

  command {
    salmon --no-version-check quant -i ~{index}  --numBootstraps ~{bootstraps} --threads ~{threads} -l A --seqBias --gcBias -o quant_~{run} \
    ~{if(is_paired) then "-1 " + reads[0] + " -2 "+ reads[1] else "-r " + reads[0]}
  }
  # --validateMappings --rangeFactorizationBins ~{rangeFactorizationBins}

  runtime {
    docker: "combinelab/salmon:0.12.0"
    maxRetries: 3
  }

  output {
    File out = "quant_" + run
    File lib = out + "/" + "lib_format_counts.json"
    File quant = out + "/" + "quant.sf"
  }
}


task copy {
    input {
        Array[File] files
        String destination
    }

    command {
        mkdir -p ~{destination}
        cp -L -R -u ~{sep=' ' files} ~{destination}
    }

    output {
        Array[File] out = files
    }
}