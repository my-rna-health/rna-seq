workflow genome_assembly {

    String name
    Array[File] reads
    Int memory
    Int threads
    Int k = 64
    Array[File] mate_1
    Array[File] mate_2


    #call test{
    #    input:
    #        name = name, reads = reads, k = 25, memory = 2, threads = threads
    #}

    call abyss2 {
        input:
            name = name, reads = reads, k = k, memory = memory, threads = threads, mate_1 = mate_1, mate_2 = mate_2
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
        /usr/local/bin/abyss-pe k=${k} name=${name} in='${sep=" " reads}' B=${memory}G H=4 kc=${kc}
    }

    #/usr/local/bin/abyss-pe k=${k} name=${name} in='${sep=" " reads}' B=${memory}G H=4 kc=${kc}


    runtime {
        docker: "bcgsc/abyss@sha256:2e0e23a2df1e9cacd2ebbb2ddda5a29d67aca7550d95d157485fbae9a0e291b1"
    }
}

task abyss {
    String name
    Array[File] reads
    Array[File] mates
    Int k
    Int memory
    Int threads
    Int kc = 3

    command {
        /usr/local/bin/abyss-pe k=${k} name=${name} in='${sep=" " reads}' B=${memory}G H=4 kc=${kc}
    }

 runtime {
    docker: "bcgsc/abyss@sha256:2e0e23a2df1e9cacd2ebbb2ddda5a29d67aca7550d95d157485fbae9a0e291b1"
 }
}

task abyss2 {
    String name
    Array[File] reads
    Array[File] mate_1
    Array[File] mate_2
    Int k
    Int memory
    Int threads
    Int kc = 3

    command {
        /usr/local/bin/abyss-pe k=${k} name=${name} in='${sep=" " reads}' \
        mp='mpc mpd' mpc='${sep=" " mate_1}' mpd='${sep=" " mate_2}' \
        B=${memory}G H=4 kc=${kc}
    }

 runtime {
    docker: "bcgsc/abyss@sha256:2e0e23a2df1e9cacd2ebbb2ddda5a29d67aca7550d95d157485fbae9a0e291b1"
 }
}
