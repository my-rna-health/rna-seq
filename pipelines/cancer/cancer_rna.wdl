workflow cancer_rna{

}

task star_fusion {
    input {

    }

    command {

    }

    runtime {
        docker: "quay.io/biocontainers/star-fusion@sha256:e0c239d18a421dd11742dbe382d2bd6688e9fb9bcafb19044aea09cbec285a3a" #1.6.0--0
    }
}