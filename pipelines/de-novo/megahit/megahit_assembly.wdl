version development

workflow megahit_assembly {
    input {
        Array[File] reads
        String folder
    }

    call fastp { input: reads = reads }

    call copy as copy_report {
         input:
            destination = folder + "/report",
            files = [fastp.report_json, fastp.report_html]
        }

    call megahit {
        input: reads = fastp.reads_cleaned
    }

    call copy as copy_assembly{
         input:
            destination = folder,
            files = [megahit.out]
        }

    call scaffolding {
        input: contigs = megahit.contigs
    }

    call copy as copy_scaffolding{
         input:
            destination = folder + "/scaffolds",
            files = scaffolding.out
        }


    output {
        File contigs = copy_assembly.out[0]
        File scaffolds = copy_scaffolding.out[0]
    }
}

task megahit {
    input {
        Array[File] reads
        Int k_min = 31
    }
    Boolean is_paired = if(length(reads) > 1) then true else false
    #String reads_str = if(is_paired == true) then " -1 " + reads[0] + " -2 " + reads[1] else " -r " + reads[0]

    command {
        megahit --no-mercy --k-min ~{k_min} -1 ~{reads[0]} -2 ~{reads[1]}
    }

    output {
        File out = "megahit_out"
        File contigs = "megahit_out/final.contigs.fa"
    }

    runtime {
     	#docker: "vout/megahit:latest"
        docker: "quay.io/biocontainers/megahit@sha256:14027ab5fe760d38debe85e002d783e72e6037dda039cba5b5272fa8aa7cd76e" #1.2.8--h8b12597_0
    }
}

task quast {

    input {
        Array[File]+ contigs
        File? reference
        Int? threads = 4
        File? features
        String output_folder = "results"
        Int min_contig = 50
        String? type
    }

    command {
        quast.py ~{if defined(reference) then "--reference " + reference else ""} \
         ~{if defined(threads) then "--threads " + threads else ""} ~{sep=" " contigs} \
         --output ~{output_folder} \
         ~{if defined(features) then "--features " + features + (if(defined(type)) then "--type " + type else "") else "" } \
         --min-contig ~{min_contig} \
         ~{sep=" " contigs}
    }

    runtime {
        docker: "quay.io/biocontainers/quast@sha256:89c337541c3bc92bed901b6215231a5b6f18bed86e25b5f94a97fee73d0e7c13" #5.0.2--py27pl526ha92aebf_0
    }

    output {
        File out = output_folder
    }
}



task scaffolding {
    input {
        File contigs
        Int k = 75
    }

    command {
        /opt/SOAPdenovo2-master/fusion/SOAPdenovo-fusion -S -D -K ~{k} -c ~{contigs} -g scaffold
    }

    output {
        File contig = "scaffold.contig"
        File arc = "scaffold.Arc"
        File contig_index = "scaffold.ContigIndex"
        File conver = "scaffold.conver"
        File preGraph= "scaffold.preGraphBasic"
        File edge = "scaffold.updated.edge"
        Array[File] out = [contig, contig_index, arc, conver, preGraph,edge]
    }

     runtime {
         	#docker: "vout/megahit:latest"
            docker: "pegi3s/soapdenovo2"
        }
}

task fastp {
    input {
        Array[File] reads
    }

    Boolean is_paired = if(length(reads) > 1) then true else false

    command {
        fastp --cut_front --cut_tail --cut_right --trim_poly_g --trim_poly_x --overrepresentation_analysis \
            -i ~{reads[0]} -o ~{basename(reads[0], ".fq.gz")}_cleaned.fq.gz \
            ~{if( is_paired ) then "--detect_adapter_for_pe " + "--correction -I "+reads[1]+" -O " + basename(reads[1], ".fq.gz") +"_cleaned.fq.gz" else ""}
    }

    runtime {
        docker: "quay.io/biocontainers/fastp@sha256:ac9027b8a8667e80cc1661899fb7e233143b6d1727d783541d6e0efffbb9594e" #0.20.0--hdbcaa40_0
    }

    output {
        File report_json = "fastp.json"
        File report_html = "fastp.html"
        Array[File] reads_cleaned = if( is_paired )
            then [basename(reads[0], ".fq.gz") + "_cleaned.fq.gz", basename(reads[1], ".fq.gz") + "_cleaned.fq.gz"]
            else [basename(reads[0], ".fq.gz") + "_cleaned.fq.gz"]
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