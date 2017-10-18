import "quality.wdl" as quality
import "quantification.wdl" as quant

workflow rna_seq {

    File conditions
    File samplesFile
    File transcriptomeIndex

    Array[Array[String]] series = read_tsv(conditions)

    scatter (condition in series) {
     call quality.cleanup as cleaning  { input:
        condition = condition
      }

      call quant.quantify as counting { input:
          transcriptome = transcriptomeIndex,
          samples = cleaning.out["trimmed"]
      }
    }

    output {
        Array[File] out = counting.out
    }

}
