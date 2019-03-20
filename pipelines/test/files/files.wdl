version development

workflow files {
    input {
        Array[String] files
        String destination
    }

    call copy{
        input: files = files,
        destination = destination
    }

    output {
       Array[String] out = copy.out
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