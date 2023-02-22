version development
import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

workflow download_samples{
    input {
        Array[String] samples = ["SAMEA3231268","SAMEA3231287"]
        String? email
        String destination
    }

    scatter (sample in samples) {
        call download {
            input: sample = sample, email = email
        }
    }

    call files.copy as copy_results {
             input:
             destination = destination, files = download.out
    }

    output {
       Array[File] results =  copy_results.out
    }
}

task download {
    input {
        String sample
        String format = "READS_FASTQ"
        String output_location = "./"
        String protocol = "FTP"
        String? email
    }

    command {
        mkdir -p reads_fastq/~{sample}
        curl -o reads_fastq/~{sample}/~{sample}.json "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=~{sample}&result=read_run&fields=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,instrument_model,library_name,library_layout,library_strategy,library_selection,read_count,experiment_title,study_title,experiment_alias,fastq_ftp,submitted_ftp,sra_ftp,sample_title&format=json&download=true&limit=0"
        curl -o reads_fastq/~{sample}/~{sample}.tsv "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=~{sample}&result=read_run&fields=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,instrument_model,library_name,library_layout,library_strategy,library_selection,read_count,experiment_title,study_title,experiment_alias,fastq_ftp,submitted_ftp,sra_ftp,sample_title&format=tsv&download=true&limit=0"
        java -jar /opt/ena/ena-file-downloader.jar --accessions=~{sample} --format=~{format} --location=~{output_location} --protocol=~{protocol} --asperaLocation=null ~{email}
        echo "~{sample} successfully downloaded!"
    }

    runtime {
        docker: "quay.io/comp-bio-aging/ena-downloader:master"
        maxRetries: 1
    }

    output {
        File json = "reads_fastq/"+sample + "/" + sample + ".json"
        File tsv = "reads_fastq/"+sample + "/" + sample + ".tsv"
        File out = "reads_fastq/"+sample
    }
}