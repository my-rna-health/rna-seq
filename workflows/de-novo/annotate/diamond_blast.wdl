workflow Diamond_Blast {

  Int threads
  File db
  File query
  String result_name
  String results_folder
  String mode
  String output_format

  call diamond_blast {
      input:
        threads = threads,
        database = db,
        name = result_name,
        query = query,
        mode = mode,
        output_format = output_format
    }

  call copy as copy_results {
    input:
        files = [diamond_blast.out],
        destination = results_folder
  }

  output {
       File out = copy_results.out[0]
  }

}

task make_examples{

    File alignment_bam
    File alignment_bai
    File fasta
    File fai
    File example_dir

    command {
        mkdir data
        cp ${alignment_bam} ${alignment_bai} ${fai} ${fasta} data

        /opt/deepvariant/bin/make_examples \
            --mode calling \
            --ref data/${basename(fasta)} \
            --reads data/${basename(alignment_bam)} \
            --examples ${example_dir}/output.examples.tfrecord \
            --regions "chr20:10,000,000-10,010,000"
    }

    runtime {
        docker: "gcr.io/deepvariant-docker/deepvariant"
    }

    output {
        File out = example_dir+"output.examples.tfrecord"
    }
}

task diamond_blast {

  Int threads
  File database
  File query
  String name
  String mode
  String output_format

    command {
        diamond ${mode} -d ${database}  -q ${query} \
          --more-sensitive -o ${name}.m8 \
          -f ${output_format}
     }
    #	qseqid means Query Seq - id
    #	qlen means Query sequence length
    #	sseqid means Subject Seq - id
    #	sallseqid means All subject Seq - id(s), separated by a ';'
    #	slen means Subject sequence length
    #	qstart means Start of alignment in query
    #	qend means End of alignment in query
    #	sstart means Start of alignment in subject
    #	send means End of alignment in subject
    #	qseq means Aligned part of query sequence
    #	sseq means Aligned part of subject sequence
    #	evalue means Expect value
    #	bitscore means Bit score
    #	score means Raw score
    #	length means Alignment length
    #	pident means Percentage of identical matches
    #	nident means Number of identical matches
    #	mismatch means Number of mismatches
    #	positive means Number of positive - scoring matches
    #	gapopen means Number of gap openings
    #	gaps means Total number of gaps
    #	ppos means Percentage of positive - scoring matches
    #	qframe means Query frame
    #	btop means Blast traceback operations(BTOP)
    #	staxids means unique Subject Taxonomy ID(s), separated by a ';' (in numerical order)
    #	stitle means Subject Title
    #	salltitles means All Subject Title(s), separated by a '<>'
    #	qcovhsp means Query Coverage Per HSP
    #	qtitle means Query title

  runtime {
    docker: "quay.io/comp-bio-aging/diamond:latest"
  }

  output {
       File out = name + ".m8"
  }

}


task copy {
    Array[File] files
    String destination

    command {
        mkdir -p ${destination}
        cp -L -R -u ${sep=' ' files} ${destination}
    }

    output {
        Array[File] out = files
    }
}