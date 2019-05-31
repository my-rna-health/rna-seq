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

workflow bs_map_fast {
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

   

    output {

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


task copy {
    input {
        Array[File] files
        String destination
    }

    command {
        mkdir -p ~{destination}
        cp -L -R -u ~{sep=' ' files} ~{destination}
        declare -a files=(~{sep=' ' files})
        for i in ~{"$"+"{files[@]}"};
          do
              value=$(basename ~{"$"}i)
              echo ~{destination}/~{"$"}value
          done
    }

    output {
        Array[File] out = read_lines(stdout())
    }
}