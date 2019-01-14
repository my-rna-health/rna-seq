version 1.0

workflow quantification {

    input {
        Array[String] samples
        Map[String, File] transcriptomes
        String samples_folder

        Int rangeFactorizationBins = 4
        Int bootstraps = 128
        Int extract_threads = 6

    }

    scatter(sample in samples) {
        call get_runs{
            input:
                gsm = sample
        }

        scatter(run in get_runs.zipped) {
            String srr = run.left
            Boolean is_paired = (run.right == "SINGLE")
            call download {input: sra = srr }
            call extract {input: sra = download.out, is_paired = is_paired, threads = extract_threads}
        }

    }



}


task get_runs{

    input {
        String gsm
    }

    String runinfo = gsm + "_runinfo"

    command {
        esearch -db sra -query ~{gsm} > ${gsm}.xml

        esearch -db sra -query ~{gsm} | efetch -format runinfo > ~{runinfo}.csv
        sed 's/,/\t/g' ~{runinfo}.csv > ~{runinfo}.tsv
        tail ~{runinfo}.tsv > ~{runinfo}_tail.tsv
    }


    output{
        File info = runinfo + ".tsv"
        File tail = runinfo + "_tail.tsv"
        Array[Array[String]] runinfo_tsv = read_tsv(tail)
        Array[Array[String]] runinfo_transposed = transpose(runinfo_tsv)

        Array[String] runs = runinfo_transposed[0]
        Array[String] layout = runinfo_transposed[15]
        Array[Pair[String, String]] zipped = zip(runs, layout)

    }

    runtime {
        docker: "quay.io/biocontainers/entrez-direct@sha256:a9b673ed5154f2b58b34d6f7e823b205fc8e45c5cd6883335ae107efec746c30" #10.2--pl526_0
        #maxRetries: 2
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

task download {
    input {
        String sra
    }

    command {
        download_sra.sh ~{sra}
    }

    #https://github.com/antonkulaga/biocontainers/tree/master/downloaders/sra

    runtime {
        docker: "quay.io/antonkulaga/download_sra:latest"
        #maxRetries: 2
    }

    output {
        File out = sra + ".sra"
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
