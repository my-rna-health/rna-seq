version 1.0

workflow chip {

    input {
        String treatment
        String control
        File reference
        String destination
        Boolean is_paired = true
        Int extract_threads
        Int aligner_threads
    }

    String treatment_result = destination + "/" + treatment

    String treatment_report = treatment_result + "/report"

    String treatment_alignment = treatment_result + "/alignment"

    #TREATMENT PREPROCESSING AND ALIGMENT

    call download as download_treatment {
            input:  sra = treatment
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
        input:  sra = control
    }

    call bwa as align_treatment {
        input:
            reads = fastp_treatment.reads_cleaned,
            reference = reference,
            name = "treatment_" + treatment,
            threads = aligner_threads
    }


    call samtools_conversion as convert_treatment{
        input: sam =  align_treatment.out
    }

    call samtools_sort as sort_treatment{
            input:
                bam = convert_treatment.out
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

    call bwa as align_control {
        input:
            reads = fastp_control.reads_cleaned,
            reference = reference,
            name = "control_" + control,
            threads = aligner_threads
    }


    call samtools_conversion as convert_control{
        input: sam = align_control.out
    }

    call samtools_sort as sort_control{
            input:
                bam = convert_control.out
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
    }

    command {
        download_sra.sh ~{sra}
    }

    #https://github.com/antonkulaga/biocontainers/tree/master/downloaders/sra

    runtime {
        docker: "quay.io/antonkulaga/download_sra@sha256:84369d824648d1294e5ad4cc3b365227f86923f7f6fb5a69e0b456a545a287dd"
        maxRetries: 2
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
        maxRetries: 3
    }

    output {
        Array[File] out = if(is_paired) then [prefix + "_1.fastq",  prefix + "_2.fastq"] else [prefix + ".fastq"]
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
        docker: "quay.io/biocontainers/fastp@sha256:56ca79fc827c1e9f48120cfa5adb654c029904d8e0b75d01d5f86fdd9b567bc5" #0.20.1--h8b12597_0
        maxRetries: 2
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
    }

    command {
        minimap2 -ax sr  -t ~{threads} -2 ~{reference} ~{sep=' ' reads} > ~{name}.sam
    }

    runtime {
        docker: "genomicpariscentre/minimap2@sha256:536d7cc40209d4fd1b700ebec3ef9137ce1d9bc0948998c28b209a39a75458fa"
        maxRetries: 2
      }

    output {
      File out = name + ".sam"
    }
}

task bwa {
    input {
        Array[File] reads
        File reference
        String name
        Int threads
    }

    String index_name = basename(reference, ".fa")

    command {
        bwa index -p ~{index_name} ~{reference}
        bwa mem -t ~{threads} ~{index_name} ~{sep=' ' reads} > ~{name}.sam
    }

    runtime {
        docker: "biocontainers/bwa@sha256:cb0cde0bf4e832697af16698369abc222606779237809a1b17fa92758c54d0be"
        maxRetries: 2
      }

    output {
      File out = name + ".sam"
    }
}

task samtools_conversion {
    input {
        File sam
    }

    String name = basename(sam, ".sam")

    command {
       samtools view -bS ~{sam} > ~{name}.bam
    }

    runtime {
        docker: "biocontainers/samtools@sha256:6644f6b3bb8893c1b10939406bb9f9cda58da368100d8c767037558142631cf3"
        maxRetries: 2
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
        docker: "biocontainers/samtools@sha256:6644f6b3bb8893c1b10939406bb9f9cda58da368100d8c767037558142631cf3"
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
            docker: "quay.io/biocontainers/bedtools@sha256:a0bb135afdec53be4b953a9a8efbc801cdb90706e6e63e11e3f60b06b8444f78" #2.23.0--he941832_1
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
        docker: "quay.io/biocontainers/macs2@sha256:7057ebb45ee9a185f132832a932aed54e2fc3d17bf638ae6c4b5ed201a6029d8" #2.1.1.20160309--py27h7eb728f_3
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