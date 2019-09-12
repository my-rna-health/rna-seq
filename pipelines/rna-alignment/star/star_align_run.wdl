version development

import "extract_run.wdl" as extractor

struct AlignedRun {
    String run
    File sorted
    File to_transcriptome
    File summary
    File log
    File progress
    File reads_per_gene
    File junctions
}


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
        run = run,
        reads = extract_run.out.cleaned_reads,
        index_dir = index_dir
    }

    call extractor.copy as copy_quant{
        input:
            destination = extract_run.out.folder,
            files = [star_align.out.sorted, star_align.out.to_transcriptome,
                star_align.out.summary, star_align.out.log, star_align.out.progress,
                star_align.out.reads_per_gene, star_align.out.junctions
       ]
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
        AlignedRun out = star_align.out

    }

}

task star_align {
    input {
        String run
        Array[File] reads
        Directory index_dir
        Float threshold  = 0.2
        Int threads = 4
        Float minOverLread = 0.2
    }


  command {
    /usr/local/bin/STAR \
        --runThreadN ~{threads} \
        --genomeDir ~{index_dir} \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM GeneCounts \
        --readFilesCommand gunzip -c \
        --outFilterMatchNminOverLread ~{threshold} \
        --outFilterScoreMinOverLread ~{minOverLread} \
        --readFilesIn ~{sep=" " reads}
  }

  runtime {
    docker: "quay.io/biocontainers/star@sha256:f9b0406354ff2e5ccfadaef6fde6367c7bcb4bdc7e67920f0f827a6ff6bf4fb5" #2.7.2b--0
  }

  output {
    AlignedRun out= object {
      run: run,
      sorted: "Aligned.sortedByCoord.out.bam",
      to_transcriptome: "Aligned.toTranscriptome.out.bam",
      summary: "Log.final.out",
      log: "Log.out",
      progress: "Log.progress.out",
      reads_per_gene: "ReadsPerGene.out.tab",
      junctions: "SJ.out.tab"
    }
  }

}
