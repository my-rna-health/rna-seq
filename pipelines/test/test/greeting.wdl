version development

workflow greeting {
    input {
        String name
        String intro = "Hello"
        Int delay = 1
    }

    call greet {
        input:
            name = name,
            intro = intro,
            delay = delay
    }

    output {
        String out = greet.out
    }

}

task greet {
    input {
        String name
        String intro
        Int delay
    }

    command {
        echo "~{intro} ~{name}!"
        sleep ~{delay}
        echo "Nice to meet you!"
    }

    output {
        String out = read_string(stdout())
    }
}

