workflow quantification {

    File batch
    File references
    Int threads
    File samples_folder
    String results_folder

    call prepare_samples {
        input: samples = batch,references = references, samples_folder = samples_folder
    }

    call copy as report_invalid {
        input: files = [prepare_samples.invalid], destination = results_folder + "/invalid"
    }

    Array[Array[String]] cached_samples = if(read_string(prepare_samples.cached)=="") then [] else read_tsv(prepare_samples.cached)

    Array[Array[String]] novel_samples = if(read_string(prepare_samples.novel)=="") then [] else read_tsv(prepare_samples.novel)

    scatter(row in novel_samples) {
            #"GSM",	"GSE",	"Species",	"Sequencer",
            #"Type", "Sex",	"Age",	"Tissue",
            #"Extracted molecule", "Strain",
            #"Comments", "salmon", "transcriptome", "gtf"

        String gsm = row[0]
        String gse = row[1]
        File salmon_index = row[12]

        call get_sample {
            input: gsm_id = gsm
        }

        File sample_row = get_sample.tsv

        Array[String] sample = read_tsv(sample_row)[0]
        Boolean is_paired = if(sample[1]=="paired") then true else false
        Array[File] reads = if(is_paired) then [sample[2], sample[3]] else [sample[2]]
        Array[File] to_copy = if(is_paired) then [sample_row, sample[2], sample[3], sample[4]] else [sample_row, sample[2], sample[4]]

        call copy as copy_sample{
            input:
                destination = samples_folder + "/" + gse + "/" + gsm + "/" + "raw",
                files = to_copy
        }

        call fastp{
            input: reads = reads, is_paired = is_paired
        }

        call copy as copy_sample_cleaned{
                input:
                    destination = samples_folder + "/" + gse + "/" + gsm + "/" + "cleaned",
                    files = fastp.reads_cleaned
            }

        call copy as copy_sample_report{
                input:
                    destination = samples_folder + "/" + gse + "/" + gsm + "/" + "report",
                    files = [fastp.report_json, fastp.report_html]
            }

        call salmon {
            input:
                index = salmon_index,
                reads = fastp.reads_cleaned,
                is_paired = is_paired,


        }

    }

    scatter(row in cached_samples) {

        call echo as echo_cached {
            input: message = row[0] + "_cached"
        }

    }

}

task fastp {

    Array[File] reads
    Boolean is_paired

    command {
        fastp --cut_by_quality5 --cut_by_quality3 --trim_poly_g --overrepresentation_analysis \
            -i ${reads[0]} -o ${basename(reads[0], ".fastq.gz")}_cleaned.fastq.gz \
            ${if( is_paired ) then "--correction -I "+reads[1]+" -O " + basename(reads[1], ".fastq.gz") +"_cleaned.fastq.gz" else ""}
    }

    runtime {
        docker: "quay.io/biocontainers/fastp@sha256:522832170d976e8fb70ccfbb8f5143d150ec2f36472316cc53276d3c50302c54"
    }

    output {
        File report_json = "fastp.json"
        File report_html = "fastp.html"
        Array[File] reads_cleaned = if( is_paired )
            then [basename(reads[0], ".fastq.gz") + "_cleaned.fastq.gz", basename(reads[1], ".fastq.gz") + "_cleaned.fastq.gz"]
            else [basename(reads[0], ".fastq.gz") + "_cleaned.fastq.gz"]
    }
}


task get_sample {

  String gsm_id

  command {
    /opt/geoparse/run.py --location ./ --filetype fastq --keep_sra false ${gsm_id}
  }

  runtime {
    docker: "quay.io/comp-bio-aging/geoparse@sha256:3ae3b327d6b01ed94dc41193c6dfc771edb68d58638dde8ae9fed909d3fa54d9"
  }

  output {
    File tsv = "output.tsv"
    File json = "output.json"
  }

}

task prepare_samples {
    File samples
    File references
    File samples_folder

    command {
        /scripts/run.sc --samples ${samples} --references ${references} --cache ${samples_folder}
    }

    runtime {
        docker: "quay.io/comp-bio-aging/prepare-samples@sha256:9aaa223ff520634bb0357500ffb90aa80315729e0870ebbc7da4a4b31c382a2c"
    }

    output {
        File invalid = "invalid.tsv"
        File novel = "novel.tsv"
        File cached = "cached.tsv"
    }
}


task salmon {
  File index
  Array[File] reads
  Boolean is_paired
  Int threads

  command {
    salmon quant -i ${index} --threads ${threads} -l A \
    -1 ${reads[0]} -o transcripts_quant \
    ${if(is_paired) then "-2 " + reads[1] else ""}
  }

  runtime {
    docker: "combinelab/salmon@sha256:8a5f0de02b0df1b2571f8200e276c09ef1dd499ca13a883577230d85d8e644c3"
  }

  output {
    File out = "transcripts_quant"
  }
}


task copy {
    Array[File] files
    String destination

    command {
        mkdir -p ${destination}
        cp -L -R -u ${sep=' ' files} ${destination}
    }

    output {
        Array[File] out = files
    }
}

task echo {
    String message

    command {
        echo ${message} >> /pipelines/test/echo.txt
    }

    output {
        String out = message
    }
}