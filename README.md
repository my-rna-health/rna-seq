RNA-Seq and other bioinformatics workflows
==========================================

This repository containers our pipelines, containers and scripts for them.
Pipelines are written in [WDL language](https://openwdl.org/) and Cromwell is used to run them. 

RNA-Seq pipeline
----------------

The pipeline is executed by Cromwell workflow execution engine. 
There is also a microservice that allows starting the pipeline directly from CKAN. 
The main part of the pipeline is quantification workflow. 
The quantification part requires only GSM ids and list of indexed transcriptomes as input and consists of 4 sub workflows: 
* quantification workflow gets the sample ids and calls quant_sample workflow for each of them. 
* quant_sample workflow extracts sample metadata by GSM id (with GEOfetch tool that we developed), gets SRR runs information and then starts quant_run workflow that does the quantification. 
* quant_run workflow:
    * with extract_run subworkflow:
      * downloads (with either prefetch or aspera connect) that data
      * extracts (with fasterq-dump)
      * does adapter/quality trimming and reporting (with fastp)
    * does quantification of samples with Salmon quantification software
    * copies (together with quality report and metadata) to the output folders.
Before running the pipeline te user must make sure that indexes for Salmon are built. 
For this purpose quant_index pipeline is used where an array of indexes is provided, for example::
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
To run with SRR as ids quant_by_runs.wdl is used instead of quantification, the rest remains the same:
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


Bs-Seq pipeline
---------------

To use Bs-Seq pipeline the user must first build the index, i.e:
```json5
{
  "bitmapper_index.index_name": "homo_sapiens",
  "bitmapper_index.genome": "/data/ensembl/100/species/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
  "bitmapper_index.destination": "/data/indexes/bitmapper/ensembl_100/"
}
```