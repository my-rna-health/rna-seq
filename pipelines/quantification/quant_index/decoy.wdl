version development

workflow decoy{
 input{
        File gtf
        File genome
        File transcriptome
        Int threads = 16
    }

    call generate{
        input:
            gtf = gtf, genome = genome, transcriptome = transcriptome, threads = 16
    }

    output {
        File out = generate.out
    }
}

task generate{

    input{
        File gtf
        File genome
        File transcriptome
        Int threads
    }

    command {
        /opt/SalmonTools/scripts/generateDecoyTranscriptome.sh -a ~{gtf} -g ~{genome} -t ~{transcriptome} -j ~{threads} -o output
    }

    runtime {
        docker: "quay.io/comp-bio-aging/salmon-tools"
    }

    output {
        File out = "output"
    }
}
