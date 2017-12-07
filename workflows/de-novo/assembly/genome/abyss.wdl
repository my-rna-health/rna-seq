workflow genome_assembly {

}

task abyss {
 command {
    abyss-pe k=64 name=ecoli lib='pea peb' mp='mpc mpd' \
        pea='pea_1.fa pea_2.fa' peb='peb_1.fa peb_2.fa' \
        mpc='mpc_1.fa mpc_2.fa' mpd='mpd_1.fa mpd_2.fa' \
        B=26G H=4 c=3 \


    abyss-pe name=hsapiens np=64 k=144 q=15 v=-v l=40 s=1000 n=10 \
    B=26G H=4 c=3 \
    S=1000-10000 N=7 mp6k_de=--mean mp6k_n=1 \
    lib=pe400 pe400=$(<pe400.in) \
    mp=mp6k mp6k=$(<mp6k+unknown.in)
 }

 runtime {
    docker: "bcgsc/abyss"
 }
}