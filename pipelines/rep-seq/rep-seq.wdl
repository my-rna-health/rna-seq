version development

workflow rep_seq{
    input {
        File fasta
        String name
        String destination
    }

    call mixcr {
        input: name = name, fasta = fasta
    }

    call copy {
        input: files = mixcr.reports, destination = destination
    }

    output {
        Array[File] reports = copy.out
    }
}

task mixcr {
    input {
        File fasta
        String name
        String memory = "10G"
    }

    command {
        set JAVA_OPTS="-Xms128m -Xmx~{memory}"
        mixcr analyze amplicon --species hsa --starting-material rna --5-end v-primers --3-end c-primers --adapters adapters-present ~{fasta} ~{name}
    }

    #java -Xmx~{memory} -jar /mixcr/mixcr-3.0.3/mixcr.jar analyze amplicon --species hsa --starting-material rna --5-end v-primers --3-end c-primers --adapters adapters-present ~{fasta} ~{name}

    output {
        Array[File] reports =  glob(name + ".*")
    }

    runtime {
        docker: "mgibio/mixcr"
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
