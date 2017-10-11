workflow check_dna_assembly
{

}

task quast {

	runtime {
	    docker: "quay.io/biocontainers/quast"
	}
}