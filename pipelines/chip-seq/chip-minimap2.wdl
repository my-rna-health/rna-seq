version 1.0

workflow chip {

    input {
        String treatment
        String control
        File reference
        String destination
        Boolean is_paired = true
        Int extract_threads
        Int aligner_threads = 12
        Int max_memory_gb = 28
        Boolean aspera_download = true
    }

    String treatment_result = destination + "/" + treatment

    String treatment_report = treatment_result + "/report"

    String treatment_alignment = treatment_result + "/alignment"


    #TREATMENT PREPROCESSING AND ALIGMENT

    call download as download_treatment {
            input:  sra = treatment, aspera_download = aspera_download
        }

    call extract as extract_treatment {
        input:  sra = download_treatment.out,  is_paired = is_paired, threads = extract_threads
    }


    call fastp as fastp_treatment {
        input: reads = extract_treatment.out, is_paired = is_paired
    }


    call copy as copy_report_treatment {
        input: files = [fastp_treatment.report_json, fastp_treatment.report_html], destination = treatment_report
    }


    call download as download_control {
        input:  sra = control, aspera_download = aspera_download
    }

    call minimap2 as minimap2_treatment{
            input: reads = fastp_treatment.reads_cleaned, reference = reference, name = "treatment_" + treatment, threads = aligner_threads, max_memory = max_memory_gb
        }


    call samtools_sort as sort_treatment{
            input:
                bam = minimap2_treatment.out
        }

    call coverage as coverage_treatment {
        input:
            bam = sort_treatment.out
    }

    call copy as copy_treatment_aligned {
        input: files = [sort_treatment.out, coverage_treatment.out], destination = treatment_alignment
    }


    #CONTROL PREPROCESSING AND ALIGMENT

    String control_result = destination + "/" + control

    String control_report = control_result + "/report"

    String control_alignment = control_result + "/alignment"


    call extract as extract_control {
        input:  sra = download_control.out,  is_paired = is_paired, threads = extract_threads
    }


    call fastp as fastp_control {
        input: reads = extract_control.out, is_paired = is_paired
    }

    call copy as copy_report_control {
        input: files = [fastp_control.report_json, fastp_control.report_html], destination = control_report
    }

    call minimap2 as minimap2_control{
        input:
            reads = fastp_control.reads_cleaned,
            reference = reference,
            name = "control_" + control,
            threads = aligner_threads,
            max_memory = max_memory_gb
    }


    call samtools_sort as sort_control{
            input:
                bam = minimap2_control.out
        }

    call coverage as coverage_control {
        input:
            bam = sort_control.out
    }

    call copy as copy_control_aligned {
        input: files = [sort_control.out, coverage_control.out], destination = control_alignment
    }

    #MAC2 PEAK CALLING

    call macs2 {
        input:
         treatment = [sort_treatment.out],
         control = [sort_control.out],
         outDir = "result",
        sampleName = treatment
    }

    call copy as copy_result {
        input: files = [macs2.excel, macs2.narrow_peaks, macs2.summits, macs2.model_r], destination = treatment_result
    }

}


task download {
    input {
        String sra
        Boolean aspera_download
    }
    #prefetch --ascp-path "/root/.aspera/connect/bin/ascp|/root/.aspera/connect/etc/asperaweb_id_dsa.openssh" --force yes -O results ~{sra}
    command {
        ~{if(aspera_download) then "download_sra_aspera.sh " else "prefetch --force yes -O results -t http "} ~{sra}
    }

    #https://github.com/antonkulaga/biocontainers/tree/master/downloaders/sra

    runtime {
        docker: "quay.io/comp-bio-aging/download_sra:sra_2.10.0"
        maxRetries: 1
    }

    output {
        File? a = "results" + "/" + sra + ".sra"
        File? b = "results" + "/" + sra + "/" + sra + ".sra"
        File out = select_first([a, b])
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
        docker: "quay.io/comp-bio-aging/download_sra:sra_2.10.0"
        maxRetries: 1
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
        fastp --cut_front --cut_tail --cut_right --trim_poly_g --trim_poly_x --overrepresentation_analysis \
        -i ~{reads[0]} -o ~{basename(reads[0], ".fastq.gz")}_cleaned.fastq.gz \
        ~{if( is_paired ) then "--detect_adapter_for_pe " + "--correction -I "+reads[1]+" -O " + basename(reads[1], ".fastq.gz") +"_cleaned.fastq.gz" else ""}
    }

    runtime {
        docker: "quay.io/biocontainers/fastp@sha256:56ca79fc827c1e9f48120cfa5adb654c029904d8e0b75d01d5f86fdd9b567bc5" #0.20.1--h8b12597_0
    }

    output {
        File report_json = "fastp.json"
        File report_html = "fastp.html"
        Array[File] reads_cleaned = if( is_paired )
                                    then [basename(reads[0], ".fastq.gz") + "_cleaned.fastq.gz", basename(reads[1], ".fastq.gz") + "_cleaned.fastq.gz"]
                                    else [basename(reads[0], ".fastq.gz") + "_cleaned.fastq.gz"]
    }
}

task minimap2 {
    input {
        Array[File] reads
        File reference
        String name
        Int threads
        Int max_memory
    }

    command {
        minimap2 -ax sr  -t ~{threads} -2 ~{reference} ~{sep=' ' reads} | samtools view -bS - > ~{name}.bam
    }

    runtime {
        docker: "quay.io/comp-bio-aging/minimap2@sha256:f5d43a4d857fc56bfa4e98df1049a7b9c8af0f1bf604580eb074953a00b455cd"
        maxRetries: 2
        docker_memory: "~{max_memory}G"
        docker_cpu: "~{threads+1}"
      }

    output {
      File out = name + ".bam"
    }
}


task samtools_sort {
    input {
        File bam
    }

    String name = basename(bam, ".bam")

    command {
       samtools sort ~{bam}  -o ~{name}_sorted.bam
    }

    runtime {
        docker: "biocontainers/samtools@sha256:da61624fda230e94867c9429ca1112e1e77c24e500b52dfc84eaf2f5820b4a2a" #v1.9-4-deb_cv1
        maxRetries: 2
      }

    output {
        File out = name + "_sorted.bam"
      }
}

task coverage {
    input {
        File bam
    }

    String name = basename(bam, ".bam")

    command {
        bedtools genomecov -bg -ibam ~{bam} > ~{name}.bedgraph
    }

     runtime {
            docker: "quay.io/biocontainers/bedtools@sha256:02e198f8f61329f9eafd1b9fc55828a31020b383403adec22079592b7d868006" #2.29.2--hc088bd4_0
            maxRetries: 2
          }

    output {
        File out = name + ".bedgraph"
    }
}


task macs2 {

    input{
        Array[File] treatment
        Array[File] control
        String outDir
        String sampleName
    }

    command {
        macs2 callpeak \
        --treatment ~{sep=' ' treatment} \
        --control ~{sep=' ' control} \
        --name ~{sampleName}
    }

    runtime {
        docker: "quay.io/biocontainers/macs2@sha256:7057ebb45ee9a185f132832a932aed54e2fc3d17bf638ae6c4b5ed201a6029d8" #2.1.1.20160309--py27h7eb728f_3
        maxRetries: 2
    }

    output {
        #see + for more info
       String suffix = "_peaks"
       File excel = sampleName + suffix + ".xls"
       File narrow_peaks = sampleName + suffix + ".narrowPeak"
       #File broad_peaks = sampleName + suffix + ".broadPeak"
       #File gapped_peaks = sampleName + suffix + ".gappedPeak"
       File summits = sampleName + "_summits.bed"
       File model_r = sampleName + "_model.r"
    }
}

task macs2_simple {

    input{
        Array[File] treatment
        String outDir
        String sampleName
    }

    command {
        macs2 callpeak \
        --broad --call-summits
        --treatment ~{sep=' ' treatment} \
        --name ~{sampleName}
    }

    runtime {
        docker: "quay.io/biocontainers/macs2@sha256:480dd8e83edf36bbd53117c944a5fc3f5e707aaac981496bfad8da3bfd269711" #2.2.7.1--py38h0213d0e_1
        maxRetries: 2
    }

    output {
        #see https://github.com/taoliu/MACS for more info
        String suffix = "_peaks"
        File excel = sampleName + suffix + ".xls"
        File narrow_peaks = sampleName + suffix + ".narrowPeak"
        #File broad_peaks = sampleName + suffix + ".broadPeak"
        #File gapped_peaks = sampleName + suffix + ".gappedPeak"
        File summits = sampleName + "_summits.bed"
        File model_r = sampleName + "_model.r"
    }
}



task copy {
    input {
        Array[File] files
        String destination
    }

    String where = sub(destination, ";", "_")

    command {
        mkdir -p ~{where}
        cp -L -R -u ~{sep=' ' files} ~{where}
        declare -a files=(~{sep=' ' files})
        for i in ~{"$"+"{files[@]}"};
        do
        value=$(basename ~{"$"}i)
        echo ~{where}/~{"$"}value
        done
    }

    output {
        Array[File] out = read_lines(stdout())
    }
}