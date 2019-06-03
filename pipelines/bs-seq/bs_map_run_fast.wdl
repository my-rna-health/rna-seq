version development

import "extract_run.wdl" as getter

struct MappedRun {
    String run
    String folder
    Boolean is_paired
    Array[File] report
    File cpg
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

    call getter.extract_run as extract_run {
        input:
            layout = layout,
            run = run,
            folder = output_folder,
            copy_cleaned = copy_cleaned,
            extract_threads = extract_threads
    }


    call methylation_search {
        input:
            index_folder = genome_index,
            reads = extract_run.out.cleaned_reads,
            is_paired = extract_run.out.is_paired,
            filename = run,
            threads = map_threads
    }

    call methylation_extraction {
        input:
             index_folder = genome_index,
             bmms = methylation_search.out,
             reads = extract_run.out.cleaned_reads,
             is_paired = extract_run.out.is_paired
    }

    output {
        MappedRun out = object
                           {
                               run: run,
                               folder: extract_run.out.folder,
                               is_paired: extract_run.out.is_paired,
                               report: extract_run.out.report,
                               cpg:  methylation_extraction.out
                           }
    }

}


task methylation_search {
   input {
        File index_folder
        Array[File] reads
        Boolean is_paired
        String filename
        Int threads
   }
   command {
        mkdir output
        /opt/BitMapperBS/bitmapperBS --search ~{index_folder} ~{if(is_paired) then " --seq1 " + reads[0] + " --seq2 "+ reads[1] + " --sensitive --pe" else " --seq1 " + reads[0]} -t ~{threads} --methy_out --bmm_folder output
   }

  runtime {
    docker: "quay.io/comp-bio-aging/bit_mapper_bs:master"
  }

  output {
    File out = "output"
  }
}

task methylation_extraction {

    input {
        File index_folder
        File bmms
        Array[File] reads
        Boolean is_paired
    }

    command {
       /opt/BitMapperBS/bitmapperBS --methy_extract ~{index_folder} --bmm_folder ~{bmms} --seq  ~{if(is_paired) then " --seq1 " + reads[0] + " --seq2 "+ reads[1] + " --sensitive --pe" else " --seq1 " + reads[0]}
    }

    runtime {
        docker: "quay.io/comp-bio-aging/bit_mapper_bs:master"
    }


    output {
        File out = "output_CpG.bedGraph"
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