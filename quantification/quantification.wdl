workflow quantification {

    File batch
    File references
    Int threads
    File samples_folder


    Boolean keep_sra = true

    call prepare_samples {
        input: samples = batch, references = references, samples_folder = samples_folder
    }

    call copy as report_invalid {
        input: files = [prepare_samples.invalid], destination = samples_folder +  "/batches/" +  basename(batch, ".tsv")
    }

    Array[Array[String]] cached_samples = if(read_string(prepare_samples.cached)=="") then [] else read_tsv(prepare_samples.cached)

    Array[Array[String]] novel_samples = if(read_string(prepare_samples.novel)=="") then [] else read_tsv(prepare_samples.novel)

    scatter(row in novel_samples) {
            #"GSM",	"GSE",	"Species",	"Sequencer",
            #"Type", "Sex",	"Age",	"Tissue",
            #"Extracted molecule", "Strain",
            #"Comments", "salmon", "transcriptome", "gtf"

        File salmon_index = row[12]

        call get_sample {
            input:
                gsm_id = row[0],
                gse_id = row[1],
                samples_folder = samples_folder,
                keep_sra = keep_sra
        }      

        call copy as copy_sample_novel{
            input:
                destination = get_sample.sample_raw,
                files = get_sample.to_copy
        }

        call fastp as fastp_novel{
            input: reads = get_sample.reads, is_paired = get_sample.is_paired
        }

        call copy as copy_sample_novel_cleaned{
                input:
                    destination = get_sample.sample_cleaned,
                    files = fastp_novel.reads_cleaned
            }

        call copy as copy_sample_novel_report{
                input:
                    destination = get_sample.sample_report,
                    files = [fastp_novel.report_json, fastp_novel.report_html]
            }

        call salmon as salmon_novel {
            input:
                index = salmon_index,
                reads = fastp_novel.reads_cleaned,
                is_paired = get_sample.is_paired,
                threads = threads
        }

         call copy as copy_novel_quant{
                        input:
                            destination = get_sample.sample_quant,
                            files = [salmon_novel.out]
                    }

    }

    call join_files as joined_files_novel {
        input:
            first = novel_samples,
            second = get_sample.sample,
            where = samples_folder + "/batches/" + basename(batch, ".tsv"),
            name = "samples.tsv"
    }

    output {
        File out = joined_files_novel.out
    }

}

task concat_files {
    Array[Array[String]] first
    Array[Array[String]] second
    String where
    String name

    command {
        mkdir -p ${where}
        cat ${first} >> ${where}/${name}
        cat ${second} >> ${where}/${name}
    }

    output {
        File out = where + "/" + name
    }
}


task join_files {
    Array[Array[String]] first
    Array[Array[String]] second
    String where
    String name

    command {
        mkdir -p ${where}
        join -t \t ${write_tsv(first)} ${write_tsv(second)} > ${where}/${name}
    }

    output {
        File out = where + "/" + name
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
  String gse_id
  File samples_folder
  Boolean keep_sra

  command {
    /opt/geoparse/run.py --location ./ --filetype fastq --keep_sra ${keep_sra} --header false ${gsm_id}
  }

  runtime {
    docker: "quay.io/comp-bio-aging/geoparse@sha256:03ec47e17832214f45ce9cc82dbdc691e280d4cfc9a17c36bfe9ec019d8c9aa9"
  }

  output {
    File tsv = "output.tsv"
    File json = "output.json"
    File sample_row = tsv
    Array[String] sample = read_tsv(sample_row)[0]
    Boolean is_paired = if(sample[1]=="paired") then true else false
    Array[File] reads = if(is_paired) then [sample[2], sample[3]] else [sample[2]]
    Array[File] to_copy = if(keep_sra) then
            if(is_paired) then [sample_row, sample[2], sample[3], sample[4]] else [sample_row, sample[2], sample[4]]
        else
            if(is_paired) then [sample_row, sample[2], sample[3]] else [sample_row, sample[2]]
    String sample_destination = samples_folder + "/" + gse_id + "/" + gsm_id
    String sample_raw  = sample_destination + "/" + "raw"
    String sample_cleaned = sample_destination + "/" + "cleaned"
    String sample_report = sample_destination + "/" + "report"
    String sample_quant = sample_destination + "/" + "quant"

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
        docker: "quay.io/comp-bio-aging/prepare-samples@sha256:14181c68503e307ad3f077012f2c13cec378771b5a43deee0fa1764e0dce45bd"
        #quay.io/comp-bio-aging/prepare-samples@sha256:42edc1440cd2b016083eb880824cdc7bcf3894d4d52de9c053d8391f81d062c3
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