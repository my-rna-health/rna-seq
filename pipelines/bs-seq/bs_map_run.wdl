version development

import "bs_extract_run.wdl" as getter

struct MappedRun {
    String run
    String folder
    Boolean is_paired
    Array[File] report
    File mapstats
    File aligned
    File cpg
    File counts
}

struct IndexedBamFile {
    File file
    File index
    File md5sum
}

workflow bs_map {
    input {
        String layout = "PAIRED"
        String run
        String output_folder
        File genome_index
        File genome
        Boolean copy_cleaned = false
        Int map_threads = 8
        Int extract_threads = 4
    }

    call getter.bs_extract_run as extract_run {
        input:
            layout = layout,
            run = run,
            folder = output_folder,
            copy_cleaned = copy_cleaned,
            extract_threads = extract_threads
    }


    call bitmapper {
        input:
            index_folder = genome_index,
            reads = extract_run.out.cleaned_reads,
            is_paired = extract_run.out.is_paired,
            filename = run,
            threads = map_threads
    }

    call picard_readgroups_sort {
        input: bam = bitmapper.out,
                    filename = run
    }

    call picard_mark_duplicates {
        input:
            bam = picard_readgroups_sort.out,
            outputBamPath = run + ".bam",
            metricsPath = run + ".metrics"
    }


    call copy as copy_bam {
        input:
            files = [picard_mark_duplicates.out.file, picard_mark_duplicates.out.index, picard_mark_duplicates.out.md5sum],
            destination = output_folder
    }

    call methyldackel {
        input:
            bam = picard_mark_duplicates.out.file,
            index = picard_mark_duplicates.out.index,
            genome = genome,
            threads = extract_threads
    }

    call copy as copy_methylation {
            input:
                files = [  methyldackel.cpg, methyldackel.counts],
                destination = output_folder
        }

    output {
        MappedRun out = object
        {
            run: run,
            folder: extract_run.out.folder,
            is_paired: extract_run.out.is_paired,
            report: extract_run.out.report,
            mapstats: bitmapper.stats,
            aligned: picard_readgroups_sort.out,
            cpg: methyldackel.cpg,
            counts: methyldackel.counts
        }
    }
    

}


task bitmapper {
   input {
        File index_folder
        Array[File] reads
        Boolean is_paired
        String filename
        Int threads
   }
   command {
        /opt/BitMapperBS/bitmapperBS --search ~{index_folder} ~{if(is_paired) then " --seq1 " + reads[0] + " --seq2 "+ reads[1] + " --sensitive --pe" else " --seq1 " + reads[0]} -t ~{threads} --mapstats --bam -o ~{filename}.bam
   }

  runtime {
    docker: "quay.io/comp-bio-aging/bit_mapper_bs:latest"
  }

  output {
    File out = "~{filename}.bam"
    File stats = "~{filename}.bam.mapstats"
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
        O=~{filename}_sorted.bam \
        RGID=4 \
        RGLB=lib1 \
        RGPL=illumina \
        RGPU=unit1 \
        RGSM=20 \
        SORT_ORDER=coordinate
    }

    runtime {
        docker: "quay.io/biocontainers/picard:sha256:f1b6e3b793e07529488217c80da78309d1c456315d0881fc57d414b6d3cac9d1"
    }

    output {
        File out = "~{filename}_sorted.bam"
    }

}


task picard_mark_duplicates {
    input {
        File bam
        String outputBamPath
        String metricsPath

        Int memory = 4
        Float memoryMultiplier = 3.0
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputBamPath})
        picard -Xmx~{memory}G \
        MarkDuplicates \
        INPUT=~{bam} \
        OUTPUT=~{outputBamPath} \
        METRICS_FILE=~{metricsPath} \
        VALIDATION_STRINGENCY=SILENT \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
        CLEAR_DT="false" \
        CREATE_INDEX=true \
        ADD_PG_TAG_TO_READS=false \
        CREATE_MD5_FILE=true
    }

    output {
        IndexedBamFile out = object {
          file: outputBamPath,
          index: sub(outputBamPath, ".bam$", ".bai"),
          md5: outputBamPath + ".md5"
        }
        File metricsFile = metricsPath
    }

    runtime {
        docker: "quay.io/biocontainers/picard:sha256:f1b6e3b793e07529488217c80da78309d1c456315d0881fc57d414b6d3cac9d1"
        memory: ceil(memory * memoryMultiplier)
    }
}

task methyldackel {
    input {
        File bam
        File index
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

task copy {
    input {
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