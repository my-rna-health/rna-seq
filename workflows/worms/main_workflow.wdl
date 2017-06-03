import "indexes.wdl" as index
import "quality.wdl" as quality
import "quantification.wdl" as quant
import "alignment.wdl" as star

workflow worms_single {

  File samplesFile
  File genomeFolder
  File transcriptome

  Array[Array[File]] samples = read_tsv(samplesFile)

  call index.indexes as salmon_index { input:
    genomeFolder = genomeFolder,
    transcriptomeFile = transcriptome,
    indexName = "worms"
  }

  call quality.cleanup as trimmed { input: samples = samples }

  call quant.quantify as counting { input:
    transcriptome = salmon_index.out,
    samples = trimmed.out["trimmed"]
  }

  output {
    Array[File] out = counting.out
  }

}
