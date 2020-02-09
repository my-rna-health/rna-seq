version development

workflow sspace {
    input {
        File contigs #contigs.fasta
        File libraries #libraries.txt
        Int extension = 1
        Int minimal_overlap = 32 #Minimum number of overlapping bases with the seed/contig during overhang consensus build up (default -m 32)
        Int o =20 #Minimum number of reads needed to call a base during an extension (default -o 20)
        Float r = 0.9 #Minimum base ratio used to accept a overhang consensus base (default -r 0.9)
    }

    call sspace{
        input:
            contigs = contigs,
            libraries = libraries,
            extension = extension,
            minimal_overlap = minimal_overlap,
            o = o,
            r = r
    }

    output {

    }

}

task sspace {
    input {
        File contigs #contigs.fasta
        File libraries #libraries.txt
        Int extension = 1
        Int minimal_overlap = 32 #Minimum number of overlapping bases with the seed/contig during overhang consensus build up (default -m 32)
        Int o =20 #Minimum number of reads needed to call a base during an extension (default -o 20)
        Float r = 0.9 #Minimum base ratio used to accept a overhang consensus base (default -r 0.9)
    }

    command {
        perl SSPACE_Basic.pl -l ~{libraries} -s ~{contigs} -x ~{extension} -m ~{minimal_overlap} -o ~{o} -t 0 -k 5 -a 0.70 -n 15 -p 0 -v 0 -z 0 -g 0 -T 1 -b standard_out
    }

     runtime {
        docker: "biocontainers/sspace:v2.1.1dfsg-4-deb_cv1"
     }
}