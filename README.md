RNA-Seq and other bioinformatics workflows
==========================================

This repository containers our pipelines, containers and scripts for them.
Pipelines are written in [WDL language](https://openwdl.org/) and Cromwell workflow execution engine is used to run them. 

Running pipelines
-----------------

There are three major ways of running any of the pipelines in the repository:
* with [CromwellClient](https://github.com/antonkulaga/cromwell-client) and Cromwell in server mode (recommended way). Documented at [CromwellClient](https://github.com/antonkulaga/cromwell-client) readme file.
* directly from Swagger API with Cromwell in a server mode
* with Cromwell or any other WDL-compartible tool in the console. Documented at [Official cromwell documentation](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/#step-3-running-the-workflow)

For RNA-Seq pipeline there is also an option to start the pipeline inside the CKAN system.

In order to run the pipeline three important components must be present: main workflow, subworkflows (if exist) and input yaml or json file.
For example, for RNA-Seq workflow the user should choose pipelines/quantification/quantification.wdl as main workflow and quant_sample.wdl, quant_run.wdl, extract_run.wdl as subworkflows.
In the input json GSM ids should be specified together with pates to output folder and salmon indexes. With CromwellClient it is possible to put default json values for each of te pipelines to reduce the number of input parameters.

All the tools used by the workflows are used in the from of official (at dockerhub) or unofficial (https://quay.io/repository/comp-bio-aging/)[https://quay.io/repository/comp-bio-aging/], build specifically for the pipeline, containers with dockerfile sources at [https://github.com/antonkulaga/biocontainers](https://github.com/antonkulaga/biocontainers).
There is no need to install any tools separately, only running docker daemon is needed. Docker installation [is documented](https://docs.docker.com/get-docker/) at the official Docker website.

RNA-Seq pipeline
----------------

The pipeline is executed by Cromwell workflow execution engine. 
There is also a microservice that allows starting the pipeline directly from CKAN. The main part of the pipeline is a quantification workflow. 
The quantification part requires only GSM ids and list of indexed transcriptomes as input and consists of 4 sub workflows: 
* quantification workflow gets the sample ids and calls quant_sample workflow for each of them. 
* quant_sample workflow extracts sample metadata by GSM id (with GEOfetch tool that we developed), gets SRR runs information and then starts quant_run workflow that does the quantification. 
* quant_run workflow:
    * with samples.wdl workflow:
    downloads metadata and find all sequence runs for each of the GSM ids.
    * with extract_run subworkflow:
      * downloads (with either prefetch or aspera connect) that data
      * extracts (with fasterq-dump)
      * does adapter/quality trimming and reporting (with fastp)
    * does quantification of samples with Salmon quantification software (quant_run.wdl)
    * copies (together with quality report and metadata) to the output folders.
Before running the pipeline the user must make sure that indexes for Salmon quantification tool are present and provided.
To build indexes quant_index pipeline is used where an array of transcriptomes and genomes is provided as input, for example:
```
{
  "quant_index_batch.indexes_folder": "/data/indexes/salmon/1.1.0/ensembl_99",
  "quant_index_batch.threads_per_index": 5,
  "quant_index_batch.references": [
    {
      "species": "acanthochromis_polyacanthus",
      "genome": "/data/ensembl/99/species/acanthochromis_polyacanthus/Acanthochromis_polyacanthus.ASM210954v1.dna.toplevel.fa",
      "transcriptome": "/data/ensembl/99/species/acanthochromis_polyacanthus/Acanthochromis_polyacanthus.ASM210954v1.cdna.all.fa",
      "version": "ASM210954v1",
      "subversion": "ensembl_99"
    }],
  "quant_index_batch.memory_per_index": "24G"
}
```
There are two ways of running the pipeline with quantification.wdl and GSM ids or with quant_by_runs.wdl and SRR idds.
To run Quantification pipeline the user has to provide an input json file with samples, title, salmon indexes and transcripts2genes files (optional)
```json
{
  "quantification.samples_folder": "/data/samples",  
  "quantification.experiments": [
    "GSM1580881",
    "GSM1580889"
  ],
  "quantification.title": "Human brain",
  "quantification.salmon_indexes" : {
     "Homo sapiens" : "/data/indexes/salmon/13.1/Homo_sapiens/gencode.v30",
     "Mus musculus" : "/data/indexes/salmon/14.1/Mus_musculus/gencode.vM21"    
  },
  "quantification.transcripts2genes" : {  
    "Homo sapiens" : "/data/ensembl/96/tx2gene/Homo_sapiens/hsapiens_96_tx2gene.tsv",
    "Mus musculus" : "/data/ensembl/96/tx2gene/Mus_musculus/mmusculus_96_tx2gene.tsv"    
  }
}
```
To run with SRR as ids quant_by_runs.wdl is used instead of quantification, also samples.wdl is not needed, the rest remains the same as with quantification pipeline:
```json5
{ 
  "quant_by_runs.samples_folder": "/data/samples",
  "quant_by_runs.runs": [
    "SRR1822398", "SRR1822406", "SRR3109726", "SRR3109724", "SRR489494"
  ],
  "quant_by_runs.salmon_indexes" : {
     "Homo sapiens" : "/data/indexes/salmon/14.1/Homo_sapiens/gencode.v30",
     "Mus musculus" : "/data/indexes/salmon/14.1/Mus_musculus/gencode.vM21"    
  },
  "quant_by_runs.transcripts2genes" : {  
    "Homo sapiens" : "/data/ensembl/100/tx2gene/Homo_sapiens/hsapiens_100_tx2gene.tsv",
    "Mus musculus" : "/data/ensembl/100/tx2gene/Mus_musculus/mmusculus_100_tx2gene.tsv"    
  }
}
```


Bisulfite sequencing (Bs-Seq) pipeline
-----------------------------

To use Bs-Seq pipeline the user must first build the index, i.e:
```json5
{
  "bitmapper_index.index_name": "homo_sapiens",
  "bitmapper_index.genome": "/data/ensembl/100/species/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
  "bitmapper_index.destination": "/data/indexes/bitmapper/ensembl_100/"
}
```
After building te index, the user can run bs_map_fast.wdl pipeline with extract_run.wdl as subpipeline (does the same job as extract.wdl in RNA-Seq)
Bs-Seq pipeline receives SRA run id with some parameters as output folder, layout, index path and threads as input, i.e:
```
{
  "bs_map_fast.run": "SRR948855",
  "bs_map_fast.output_folder": "/data/samples/bs-seq/SRR948855",
  "bs_map_fast.layout": "PAIRED",
  "bs_map_fast.genome_index": "/data/indexes/bitmapper/ensembl_97/homo_sapiens",
  "bs_map_fast.genome": "/data/ensembl/97/species/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
  "bs_map_fast.extract_threads": 4,
  "bs_map_fast.copy_cleaned": false,
  "bs_map_fast.map_threads": 8
}
```
It does quality control, reads deduplication and alignment as well as methylation extraction.

Other pipelines
---------------
The repository also contains functional Chip-Seq ( pipelines/chip-seq ) and de-novo DNA assembly ( pipelines/de-novo/dna) as well as tools section with WDL scripts for multiple tools. 
DNA-Seq pipeline is located in a separate repository at [github.com/antonkulaga/dna-seq](http://github.com/antonkulaga/dna-seq)
