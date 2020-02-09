version development

workflow greeting {
    input {
        String name
    }

    call greet {
        input:
            name = name
    }

    output {
        String out = greet.out
    }

}

task greet {
    input {
        String name
    }

    command {
        echo "Hello ~{name}!"
    }

    output {
        String out = read_string(stdout())
    }
}

