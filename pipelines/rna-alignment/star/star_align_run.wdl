version development

import "extract_run.wdl" as extractor

workflow star_align_run {
    input {
        String run
        String layout
        Directory index_dir
        String folder
        File tx2gene
        Map[String, String] metadata = {"run": run, "layout": layout}

        String key = "0a1d74f32382b8a154acacc3a024bdce3709"
        Int extract_threads = 4
        Int salmon_threads = 4
        Int bootstraps = 128
        Boolean copy_cleaned = false
        String prefix = ""
        Boolean aspera_download = true
    }

    call extractor.extract_run as extract_run{
        input:
            layout = layout,
            run =  run,
            folder = folder,
            copy_cleaned = copy_cleaned,
            extract_threads = extract_threads,
            aspera_download = aspera_download
    }


  call star_align {
      input:
        reads = extract_run.out.cleaned_reads,
        index_dir = index_dir
    }

    call extractor.copy as copy_quant{
    input:
       destination = extract_run.out.folder,
       files = [star_align.log, star_align.out, star_align.junctions]
    }


#    call salmon {
#        input:
#            index = salmon_index,
#            reads = extract_run.out.cleaned_reads,
#            is_paired = extract_run.out.is_paired,
#            threads = salmon_threads,
#            bootstraps = bootstraps,
#            name = prefix + run
#    }

#    call tximport {
#        input:
#            tx2gene =  tx2gene,
#            samples = salmon.quant,
#            name =  prefix + run
#    }

    output {

    }

}

task star_align {
    input {
        Array[File] reads
        Directory index_dir
        Float threshold  = 0.66
        Int threads = 4
    }


  command {
    /usr/local/bin/STAR \
        --runThreadN ~{threads} \
        --genomeDir ~{index_dir} \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM GeneCounts \
        --readFilesCommand gunzip -c \
        --outFilterMatchNminOverLread ~{threshold} \
        --readFilesIn ~{sep=" " reads} \
        --outSAMtype BAM SortedByCoordinate
  }

  runtime {
    docker: "quay.io/biocontainers/star@sha256:f9b0406354ff2e5ccfadaef6fde6367c7bcb4bdc7e67920f0f827a6ff6bf4fb5" #2.7.2b--0
  }

  output {
    File out = "Aligned.sortedByCoord.out.bam" #"Aligned.out.sam"
    File log = "Log.final.out"
    File junctions = "SJ.out.tab"
  }

}
