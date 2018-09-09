workflow quantification {

    File batch
    File references
    Int threads
    String samples_folder

    Boolean keep_sra = true
    Boolean copy_samples = false
    Boolean copy_cleaned = false
    Int rangeFactorizationBins = 4
    Int bootstraps = 10

    call prepare_samples {
        input: samples = batch, references = references, samples_folder = samples_folder
    }

    call copy as copy_prepared_samples {
        input: files = [prepare_samples.invalid, prepare_samples.novel], destination = samples_folder +  "/batches/" +  basename(batch, ".tsv")
    }

    Array[Array[String]] novel_samples = if(read_string(prepare_samples.novel)=="") then [] else read_tsv(prepare_samples.novel)

    scatter(row in novel_samples) {
            #"GSM",	"GSE",	"Species",	"Sequencer",
            #"Type", "Sex",	"Age",	"Tissue",
            #"Extracted molecule", "Strain",
            #"Comments", "salmon", "transcriptome", "gtf"

        call get_sample {
            input:
                gsm_id = row[0],
                gse_id = row[1],
                keep_sra = keep_sra
        }

        if(copy_samples) {
            call copy as copy_sample_novel{
                    input:
                        destination = samples_folder + "/" + row[1] + "/" + row[0] + "/" + "raw",
                        files = get_sample.files
                }
        }


        call process_sample {
             input:
                 sample_raw_tsv = get_sample.tsv,
                 samples_folder = samples_folder,
                 gsm_id = row[0],
                 gse_id = row[1],
                 keep_sra = keep_sra,
                 is_paired = get_sample.is_paired,
                 input_row = row
        }

        call copy as copy_sample_tsv{
            input:
                destination = process_sample.sample_destination,
                files = [process_sample.out]
        }

        call fastp as fastp_novel{
            input: reads = get_sample.files, is_paired = get_sample.is_paired
        }

        if(copy_cleaned) {
            call copy as copy_sample_novel_cleaned{
                            input:
                                destination = process_sample.sample_cleaned,
                                files = fastp_novel.reads_cleaned
                        }
        }

        call copy as copy_sample_novel_report{
                input:
                    destination = process_sample.sample_report,
                    files = [fastp_novel.report_json, fastp_novel.report_html]
            }

        call salmon as salmon_novel {
            input:
                index = row[11],
                reads = fastp_novel.reads_cleaned,
                is_paired = get_sample.is_paired,
                threads = threads,
                rangeFactorizationBins = rangeFactorizationBins,
                bootstraps = bootstraps
        }

         call copy as copy_novel_quant{
                        input:
                            destination = process_sample.sample_destination,
                            files = [salmon_novel.out]
                    }

    }

    call summarize_files {
        input:
            novel = process_sample.out,
            novel_quants = copy_novel_quant.out,
    }


    call copy as copy_concatenations {
        input: files = [summarize_files.novel_tsv], destination = samples_folder +  "/batches/" +  basename(batch, ".tsv")
    }

    output {
        File novel_tsv = summarize_files.novel_tsv
        File expressions_tsv = summarize_files.expressions_tsv
    }

}

task prepare_ml {
    Array[Array[File]] copied


    command {
        /scripts/tsv.sc concat
    }

    runtime {
        docker: "quay.io/comp-bio-aging/prepare-samples:latest"
    }

}

task summarize_files {

    Array[File] novel
    Array[Array[File]] novel_quants #just for the running order

    #TODO: find a way for docker run -v /pipelines:/pipelines quay.io/comp-bio-aging/prepare-samples process.sc update_from_json_column $(pwd)/novel.tsv $(pwd)/novel.tsv 24 25=expected_format 26=compatible_fragment_ratio


    command {
        /scripts/tsv.sc concat novel.tsv ${sep=' ' novel}
    }

    runtime {
        docker: "quay.io/comp-bio-aging/prepare-samples:latest"
    }

    output {
        File novel_tsv = "novel.tsv"
        File expressions_tsv = "novel.tsv"
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
        docker: "quay.io/biocontainers/fastp@sha256:1ae5d7ce7801391d9ed8622d7208fd7b0318a3e0c1431a039d3498d483742949" #:0.19.3--hd28b015_0
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
  Boolean keep_sra

  command {
    /opt/geoparse/run.py --location ./ --filetype fastq --keep_sra ${keep_sra} --header false ${gsm_id}
  }

  runtime {
    docker: "quay.io/comp-bio-aging/geoparse:latest"
  }

  output {
    File tsv = "output.tsv"
    File json = "output.json"

    Array[String] sample = read_tsv(tsv)[0]
    Boolean is_paired = if(sample[1]=="paired") then true else false

    Array[File] files = if(keep_sra) then
            if(is_paired) then [sample[2], sample[3], sample[4]] else [sample[2], sample[4]]
        else
            if(is_paired) then [sample[2], sample[3]] else [sample[2]]

  }

}

task process_sample {

    File sample_raw_tsv
    String samples_folder

    String gsm_id
    String gse_id
    Boolean keep_sra
    Boolean is_paired

    Array[String] input_row

    command {
        /scripts/process.sc write_sample ${write_tsv([input_row])} ${sample_raw_tsv} sample.tsv ${samples_folder} ${gse_id} false
    }

    runtime {
        docker: "quay.io/comp-bio-aging/prepare-samples:latest"
    }

    output {
        File out = "sample.tsv"
        String sample_destination = samples_folder + "/" + gse_id + "/" + gsm_id
        String sample_raw  = sample_destination + "/" + "raw"
        String sample_cleaned = sample_destination + "/" + "cleaned"
        String sample_report = sample_destination + "/" + "report"
    }

}

task prepare_samples {
    File samples
    File references
    String samples_folder

    command {
        /scripts/process.sc process --samples ${samples} --references ${references} --cache ${samples_folder}
    }

    runtime {
        docker: "quay.io/comp-bio-aging/prepare-samples:latest"
    }

    output {
        File invalid = "invalid.tsv"
        File novel = "novel.tsv"
    }
}


task salmon {
  File index
  Array[File] reads
  Boolean is_paired
  Int threads
  Int rangeFactorizationBins = 4
  Int bootstraps = 10

  command {
    salmon --no-version-check quant -i ${index}  --numBootstraps ${bootstraps} --threads ${threads} -l A --seqBias --gcBias --validateMappings --rangeFactorizationBins ${rangeFactorizationBins} -o transcripts_quant \
    ${if(is_paired) then "-1 " + reads[0] + " -2 "+ reads[1] else "-r " + reads[0]}
  }

  runtime {
    docker: "combinelab/salmon@sha256:bb9b64804d9ac79c98cc19c11a61e65bb290446beec377d46229c2686990c311" #0.11.2
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