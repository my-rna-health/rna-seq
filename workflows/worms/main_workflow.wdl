import "quality.wdl" as quality
import "quantification.wdl" as quant

workflow worms_single {

  File samplesFile
  File transcriptomeIndex

  Array[Array[File]] samples = read_tsv(samplesFile)

  call quality.cleanup as trimmed { input: samples = samples }

  call quant.quantify as counting { input:
    transcriptome = transcriptomeIndex,
    samples = trimmed.out["trimmed"]
  }

  output {
    Array[File] out = counting.out
  }

}
