workflow assemble_de_novo {


}

task assemble {

    runtime {
        docker: trinityrnaseq/trinityrnaseq

    }

}