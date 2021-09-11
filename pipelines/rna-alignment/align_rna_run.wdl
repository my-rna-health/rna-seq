version development

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files
import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/download/download_runs.wdl" as downloader
import "https://raw.githubusercontent.com/antonkulaga/rna-seq/master/pipelines/rna-alignment/star/star.wdl" as star

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
        String title = ""
        Array[String] runs
        String experiment_folder
        File index_dir
        String key = "0a1d74f32382b8a154acacc3a024bdce3709"
        Int extract_threads = 12
        Int max_memory_gb = 42
        Int align_threads = 12
        Float minOverLread = 0.2
        Float starThreshold = 0.2
        Boolean copy_extracted = true
        Boolean copy_cleaned = true
        Boolean aspera_download = true
        Boolean skip_technical = true
        Boolean original_names = false
        Boolean deep_folder_structure = true
    }

    call downloader.download_runs as download_runs{
        input:
            title = title,
            runs = runs,
            experiment_folder = experiment_folder,
            key = key,
            extract_threads = extract_threads,
            copy_cleaned = copy_cleaned,
            aspera_download = aspera_download,
            skip_technical = skip_technical,
            original_names = original_names,
            copy_extracted = copy_extracted,
    }

    Array[CleanedRun] cleaned_runs =  download_runs.out

    scatter(run in cleaned_runs) {
        String name = run.run

        call star.star_align as star_align{
            input:
                run = name,
                reads = run.cleaned_reads,
                index_dir = index_dir,
                max_memory = max_memory_gb,
                threads = align_threads,
                minOverLread = minOverLread,
                threshold = starThreshold
        }

        AlignedRun aligned =  star_align.out
        String copy_to = run.folder + "/" + "aligned"
        Pair[AlignedRun, String] pairs = (aligned, copy_to)

        call files.copy as copy_star {
            input:
            files = [aligned.sorted, aligned.to_transcriptome, aligned.summary, aligned.log, aligned.progress, aligned.reads_per_gene, aligned.junctions],
            destination =  run.folder + "/" + "aligned",
        }
    }
    output {
        Array[Pair[AlignedRun, String]] mappings = pairs
        Array[Array[File]] out = copy_star.out
    }

}