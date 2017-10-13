workflow check_dna_assembly
{
    Map[String, File] genomes

}

task quast {

    command {

    }

	runtime {
	    docker: "quay.io/biocontainers/quast:4.5--boost1.61_1"
	}

	output {
	
	}
}