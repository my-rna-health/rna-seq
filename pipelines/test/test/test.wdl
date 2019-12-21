version development

import "greeting.wdl" as greeter

workflow test {
    input {
        Array[String] experiments
        String title
        Int delay = 15
    }

     scatter(experiment in experiments) {
        call greeter.greet as greet {
                input:
                    name = experiment,
                    intro = title,
                    delay = delay
            }
     }

    output {
        Array[String] out = greet.out
    }

}