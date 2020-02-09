version development

import "greeting.wdl" as greeter

workflow hello {
    input {
        String name
    }

    call greeter.greet as greet {
        input:
            name = name
    }

    output {
        String out = greet.out
    }

}