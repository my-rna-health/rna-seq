version development

import "bs_extract_run.wdl" as getter

workflow bs_map {
    input {
        String layout = "PAIRED"
        String run
        String folder
        File genome
        Boolean copy_cleaned = false
        Int extract_threads = 4
        Int map_threads = 8
    }

    call getter.bs_extract_run as extract_run {
        input:
            layout = layout,
            run = run,
            folder = folder,
            copy_cleaned = copy_cleaned,
            extract_threads = extract_threads
    }

    call bitmapper {
        input:
            genome = genome,
            reads = extract_run.out.cleaned_reads,
            is_paired = extract_run.out.is_paired,
            filename = run,
            threads = map_threads
    }

    call picard_readgroups_sort {
        input: bam = bitmapper.out,
                    filename = run
    }


}

task bitmapper {
   input {
        File genome
        Array[File] reads
        Boolean is_paired
        String filename
        Int threads = 8
   }

   command {
        /opt/BitMapperBS/bitmapperBS --search ~{genome} ~{if(is_paired) then " --seq1 " + reads[0] + " --seq2 "+ reads[1] + " --sensitive --pe" else " --seq1 " + reads[0]} -t ~{threads} --mapstats --bam -o ~{filename}.bam
   }

  runtime {
    docker: "quay.io/comp-bio-aging/bit_mapper_bs:latest"
  }

  output {
    File out = "~{filename}.bam"
    File stats = "~{filename}.mapstats"
  }
}

task picard_readgroups_sort{
    input {
        File bam
        String filename
    }
    command {
        picard AddOrReplaceReadGroups \
        I=~{bam} \
        O=~{filename}.bam \
        RGID=4 \
        RGLB=lib1 \
        RGPL=illumina \
        RGPU=unit1 \
        RGSM=20 \
        SORT_ORDER=coordinate
    }

    runtime {
        docker: "biocontainers/picard:v2.3.0_cv3"
    }

    output {
        File out = "~{filename}.bam"
    }

}

task methyldackel {
    input {
        File bam
        File genome
        Int threads = 4
    }

    command {
        MethylDackel extract --CHH --CHG --counts -@ ~{threads} ~{genome} ~{bam}
    }


    runtime {
        docker: "quay.io/biocontainers/methyldackel@sha256:d434c3e320a40648a3c74e268f410c57649ab208fcde4da93677243b22900c55" #0.3.0--h84994c4_3
    }

    output {
        File cpg = "alignments_CpG.bedGraph"
        File counts = "alignments.counts.bedGraph"
        #File chh = "-CHH and --CHG
    }
}