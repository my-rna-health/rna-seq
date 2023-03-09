version development
import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

struct RunInfo{
    File run_folder
    File sample_folder
    Array[File] reads
    String study_accession
    String secondary_study_accession
    String sample_accession
    String secondary_sample_accession
    String experiment_accession
    String run_accession
    String tax_id
    String scientific_name
    String instrument_model
    String library_name
    String library_layout
    String library_strategy
    String library_selection
    String read_count
    String experiment_title
    String study_title
    String experiment_alias
    String fastq_ftp
    String submitted_ftp
    String sra_ftp
    String sample_title
}

workflow download_sample{
    input {
        String sample # = "SAMEA3231268"
        String? email
        String destination
    }

    call download as download_sample{
        input: sample = sample, email = email
    }

    call files.copy as copy_sample {
        input:
            destination = destination, files = [download_sample.out]
    }

    File sample_folder = copy_sample.out[0]
    scatter (run_row in download_sample.tsv_body){
        Boolean is_paired = run_row[10] == "PAIRED"
        String run = run_row[5]
        File run_folder = sample_folder + "/" + run
        Array[File] reads = if is_paired then select_all([run_folder + "/" + run + "_1.fastq.gz", run_folder + "/" + run + "_2.fastq.gz"]) else [run_folder + "/" + run + ".fastq.gz"]
        RunInfo run_info = object {
            sample_folder: sample_folder,
            run_folder: run_folder,
            reads: reads,
            study_accession: run_row[0],
            secondary_study_accession: run_row[1],
            sample_accession: run_row[2],
            secondary_sample_accession: run_row[3],
            experiment_accession: run_row[4],
            run_accession: run_row[5],
            tax_id: run_row[6],
            scientific_name: run_row[7],
            instrument_model: run_row[8],
            library_name: run_row[9],
            library_layout: run_row[10],
            library_strategy: run_row[11],
            library_selection: run_row[12],
            read_count: run_row[13],
            experiment_title: run_row[14],
            study_title: run_row[15],
            experiment_alias: run_row[16],
            fastq_ftp: run_row[17],
            submitted_ftp: run_row[18],
            sra_ftp: run_row[19],
            sample_title: run_row[20]
        }
    }

    output {
       File result_folder =  sample_folder
       Array[RunInfo] runs = run_info
    }
}

task debug {
    input {
        Array[String] to_print
    }

    command {
        echo "~{sep=" " to_print}"
    }

    output {
        Array[String] out = read_lines(stdout())
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
        head -n 1 reads_fastq/~{sample}/~{sample}.tsv > reads_fastq/~{sample}/~{sample}_header.tsv
        tail -n +2 reads_fastq/~{sample}/~{sample}.tsv > reads_fastq/~{sample}/~{sample}_body.tsv
        java -jar /opt/ena/ena-file-downloader.jar --accessions=~{sample} --format=~{format} --location=~{output_location} --protocol=~{protocol} --asperaLocation=null ~{email}
        echo "~{sample} successfully downloaded!"
    }

    runtime {
        docker: "quay.io/comp-bio-aging/ena-downloader:master"
        maxRetries: 1
    }

    output {
        File json = "reads_fastq/" + sample + "/" + sample + ".json"
        Array[Array[String]] tsv = read_tsv("reads_fastq/"+sample + "/" + sample + ".tsv")
        Array[Array[String]] tsv_header = read_tsv("reads_fastq/"+sample + "/" + sample + "_header" + ".tsv")
        Array[Array[String]] tsv_body = read_tsv("reads_fastq/"+sample + "/" + sample + "_body" + ".tsv")
        File out = "reads_fastq/"+sample
    }
}