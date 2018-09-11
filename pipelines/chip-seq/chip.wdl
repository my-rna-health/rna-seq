version 1.0

workflow chip {

    input {
        String treatment
        String control
        File reference
        String destination
        Int download_threads
        Boolean is_paired = true
    }

    String treatment_result = destination + "/" + treatment

    String treatment_report = treatment_result + "/report"

    String control_result = destination + "/" + control

    String control_report = control_result + "/report"


    call download as download_treatment {
        input:  sra = treatment, threads = download_threads, is_paired = is_paired
    }


    call fastp as fastp_treatment {
        input: reads = download_treatment.out, is_paired = is_paired
    }


    call copy as copy_report_treatment {
        input: files = [fastp_treatment.report_json, fastp_treatment.report_html], destination = treatment_report
    }


    call download as download_control {
        input:  sra = control, threads = download_threads, is_paired = is_paired
    }

    call fastp as fastp_control {
        input: reads = download_control.out, is_paired = is_paired
    }

   call copy as copy_report_control {
        input: files = [fastp_control.report_json, fastp_control.report_html], destination = control_report
    }


    call minimap2 as minimap2_treatment{
        input:
            reads = fastp_treatment.reads_cleaned,
            reference = reference
    }



    call minimap2 as minimap2_control{
        input:
            reads = fastp_control.reads_cleaned,
            reference = reference
    }

    call macs2 {
        input:
         treatment = [minimap2_treatment.out],
         control = [minimap2_control.out],
         outDir = "result",
        sampleName = treatment
    }

    call copy as copy_result {
        input: files = [macs2.excel, macs2.narrow_peaks, macs2.broad_peaks, macs2.summits, macs2.model_r], destination = treatment_result
    }

}


task download {
    input {
     String sra
     Int threads
     Boolean is_paired
    }


    # read the following explanations for parameters
    # https://edwards.sdsu.edu/research/fastq-dump/

    command {
        parallel-fastq-dump --sra-id ~{sra} --threads ~{threads} --split-files --gzip --skip-technical --readids --read-filter pass --dumpbase --clip
    }

    runtime {
        #docker: "quay.io/biocontainers/sra-tools@sha256:6d8c1daefedf6a4d00b3351d030228ca7cc4669f73f5fb76d76d1b80e79905f1"
        docker: "quay.io/biocontainers/parallel-fastq-dump@sha256:affe192bcdd2b9685ea5544200c650e22068f0acdb4249fad27b30f3673058bc" #0.6.3--py36_1
    }

    output {
        Array[File] out = if(is_paired) then [sra + "_pass_1.fastq.gz",  sra + "_pass_2.fastq.gz"] else [sra + "_pass_1.fastq.gz"]
     }
}

task fastp {

    input {
        Array[File] reads
        Boolean is_paired = true
    }

    command {
        fastp --cut_by_quality5 --cut_by_quality3 --trim_poly_g \
            -i ~{reads[0]} -o ~{basename(reads[0], ".fastq.gz")}_cleaned.fastq.gz \
            ~{if( is_paired ) then "--correction -I "+reads[1]+" -O " + basename(reads[1], ".fastq.gz") +"_cleaned.fastq.gz" else ""}
    }

    runtime {
        docker: "quay.io/biocontainers/fastp@sha256:27a12f74e8e211bdf948521060a0f6dfc44b011a2a13cd1a74daf2b126d5ce31" #:0.19.4--hd28b015_0
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
    }

    command {
        minimap2 -ax sr ~{reference} ~{sep=' ' reads} > aln.sam
    }

    runtime {
        docker: "genomicpariscentre/minimap2@sha256:536d7cc40209d4fd1b700ebec3ef9137ce1d9bc0948998c28b209a39a75458fa"
      }
    output {
      File out = "aln.sam"
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
    }

    output {
        #see https://github.com/taoliu/MACS for more info
       String suffix = "_peaks"
       File excel = sampleName + suffix + ".xls"
       File narrow_peaks = sampleName + suffix + ".narrowPeak"
       File broad_peaks = sampleName + suffix + ".broadPeak"
       File gapped_peaks = sampleName + suffix + ".gappedPeak"
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
        docker: "quay.io/biocontainers/macs2@sha256:7057ebb45ee9a185f132832a932aed54e2fc3d17bf638ae6c4b5ed201a6029d8" #2.1.1.20160309--py27h7eb728f_3
    }

    output {
        #see https://github.com/taoliu/MACS for more info
        String suffix = "_peaks"
        File excel = sampleName + suffix + ".xls"
        File narrow_peaks = sampleName + suffix + ".narrowPeak"
        File broad_peaks = sampleName + suffix + ".broadPeak"
        #File gapped_peaks = sampleName + suffix + ".gappedPeak"
        File summits = sampleName + "_summits.bed"
        File model_r = sampleName + "_model.r"
    }
}

task copy {
    input{
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