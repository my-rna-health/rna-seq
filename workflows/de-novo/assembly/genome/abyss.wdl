workflow genome_assembly {

    String name
    Array[File] reads
    Int memory
    Int threads
    Int k = 64


    call test{
        input:
            name = name, reads = reads, k = 25, memory = 2, threads = threads
    }

}

task test {
    String name
    Array[File] reads
    Int k
    Int memory
    Int threads
    Int kc = 3

    command {
        /usr/local/bin/abyss-pe k=${k} name=${name} in='${sep=" " reads}'
    }

    #/usr/local/bin/abyss-pe k=${k} name=${name} in='${sep=" " reads}' B=${memory}G H=4 kc=${kc}


    runtime {
        docker: "bcgsc/abyss@sha256:c6aed1641f2b6388a9e2864b57b7e8a0db539bfaaa9e141c28dcb2abb3b57a6c"
    }
}

task abyss {


# abyss-pe k=64 name=ecoli lib='pea peb' mp='mpc mpd' \
#         pea='pea_1.fa pea_2.fa' peb='peb_1.fa peb_2.fa' \
#         mpc='mpc_1.fa mpc_2.fa' mpd='mpd_1.fa mpd_2.fa' \
#         B=26G H=4 c=3 \


    String name
    Array[File] reads
    Array[File] mates
    Int k
    Int memory
    Int threads
    Int kc = 3

    command {
        abyss-pe name=${name} k=${k} in='${sep=" " reads}' \
          B=${memory}G H=4 kc=${kc} np=${threads}

    abyss-pe name=hsapiens np=64 k=144 q=15 v=-v l=40 s=1000 n=10 \
    B=26G H=4 c=3 \
    S=1000-10000 N=7 mp6k_de=--mean mp6k_n=1 \
    lib=pe400 pe400=$(<pe400.in) \
    mp=mp6k mp6k=$(<mp6k+unknown.in)
 }

 runtime {
    docker: "bcgsc/abyss:latest"
 }
}