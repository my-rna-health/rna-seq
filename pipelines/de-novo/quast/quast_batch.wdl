version development

import "quast.wdl" as quast


struct Quast_Parameters{
     Array[File]+ contigs
     File? reference
     String destination
     String name

}

workflow quast_batch{
    input {
       Array[Quast_Parameters] parameters
       Int? threads
    }


    scatter(param in parameters) {
        call quast.quast{
            input:
                contigs = param.contigs,
                reference = param.reference,
                destination = param.destination,
                output_name = param.name
        }
    }

}