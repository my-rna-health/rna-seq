version development

import "extract_run.wdl" as getter

struct MappedRun {
    String run
    String folder
    Boolean is_paired
    Array[File] report
    File mapstats
    File aligned
    File cpg
    File chg
    File chh
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

    call getter.extract_run as extract_run {
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


    call getter.copy as copy_bam {
        input:
            files = [picard_mark_duplicates.out.file, picard_mark_duplicates.out.index, picard_mark_duplicates.out.md5sum],
            destination = output_folder
    }

    call methyldackel {
        input:
            run = run,
            bam = picard_mark_duplicates.out.file,
            index = picard_mark_duplicates.out.index,
            genome = genome
    }

    call getter.copy as copy_methylation {
            input:
                files = [methyldackel.cpg, methyldackel.chg, methyldackel.chh],
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
            chg: methyldackel.chg,
            chh: methyldackel.chh
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
        /opt/BitMapperBS/bitmapperBS --search ~{index_folder} ~{if(is_paired) then " --seq1 " + reads[0] + " --seq2 "+ reads[1] + " --sensitive --pe" else " --seq1 " + reads[0]} -t ~{threads} --mapstats ~{filename}.mapstats --bam -o ~{filename}.bam
   }

  runtime {
    docker: "quay.io/comp-bio-aging/bit_mapper_bs:master"
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
        O=~{filename}_sorted.bam \
        RGID=4 \
        RGLB=lib1 \
        RGPL=illumina \
        RGPU=unit1 \
        RGSM=20 \
        SORT_ORDER=coordinate
    }

    runtime {
        docker: "quay.io/biocontainers/picard@sha256:93852f987fa73839e6f6eb1ecb7cfdfeb06135b789589be1a60c3dffcaf67f56" #2.20.2--0
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
          md5sum: outputBamPath + ".md5"
        }
        File metricsFile = metricsPath
    }

    runtime {
        docker: "quay.io/biocontainers/picard@sha256:93852f987fa73839e6f6eb1ecb7cfdfeb06135b789589be1a60c3dffcaf67f56" #2.20.2--0
        memory: ceil(memory * memoryMultiplier)
    }
}

task methyldackel {
    input {
        String run
        File bam
        File index
        File genome
        Boolean cytosine_report = true
        Boolean CHH = false
        Boolean CHG = false
    }

    String mode = if(cytosine_report) then "--cytosine_report" else "--counts"

    command {
        MethylDackel extract ~{mode} ~{if(CHH)  then "--CHH" else ""} ~{if(CHG)  then "--CHH" else ""}  ~{genome} ~{bam} -o $(pwd)/~{run}
    }


    runtime {
        docker: "quay.io/biocontainers/methyldackel@sha256:579532ddf7ec19a5210854a95a949d35dbde232b0a0c10129cc36b6ccea1558a" #0.4.0--hc0aa232_0
    }

    output {
        File? chg = run + "_CHG.counts.bedGraph"
        File? chh = run + "_CHH.counts.bedGraph"
        File? cpg = run + "_CpG.counts.bedGraph"
        File? report = run + ".cytosine_report.txt"
    }
}