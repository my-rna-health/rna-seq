version development

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

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

workflow star {
    input {
        String run
        Array[File] reads
        File index_dir
        Float threshold  = 0.2
        Int threads = 4
        Float minOverLread = 0.2
        Int max_memory
        String zip = "gunzip"
        String destination

    }

    call star_align{
        input: run = run, reads = reads, index_dir = index_dir, threshold = threshold, threads = threads, minOverLread = 0.2, max_memory = max_memory,
        read_files_command = if(zip == "gunzip") then "gunzip -c" else if(zip=="") then None else "bzip2 -d"
    }

    call files.copy as copy_aligned{
        input:
            destination = destination,
            files = [star_align.out.sorted, star_align.out.to_transcriptome,
                    star_align.out.summary, star_align.out.log, star_align.out.progress,
                    star_align.out.reads_per_gene, star_align.out.junctions
                    ]
    }

    output {
        AlignedRun out= object {
                            run: run,
                            sorted: copy_aligned.out[0],
                            to_transcriptome: copy_aligned.out[1],
                            summary: copy_aligned.out[2],
                            log: copy_aligned.out[3],
                            progress: copy_aligned.out[4],
                            reads_per_gene: copy_aligned.out[5],
                            junctions: copy_aligned.out[6]
                        }
    }
}

task star_align {
    input {
        String run
        Array[File] reads
        File index_dir
        Float threshold  = 0.2
        Int threads = 4
        Float minOverLread = 0.2
        Int max_memory
        String? read_files_command = "gunzip -c"
    }


    command {
        /usr/local/bin/STAR \
        --runThreadN ~{threads} \
        --genomeDir ~{index_dir} \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM GeneCounts \
        ~{"--readFilesCommand " + read_files_command} \
        --outFilterMatchNminOverLread ~{threshold} \
        --outFilterScoreMinOverLread ~{minOverLread} \
        --readFilesIn ~{sep=" " reads}
    }

    runtime {
        docker_memory: "~{max_memory}G"
        docker: "quay.io/biocontainers/star:2.7.9a--h9ee0642_0"
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