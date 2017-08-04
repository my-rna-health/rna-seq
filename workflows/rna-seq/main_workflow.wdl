import "quality.wdl" as quality
import "quantification.wdl" as quant

workflow rna_seq {

    File conditions
    File transcriptomeIndex

    Array[Array[String]] series =  read_tsv(conditions)

    scatter (condition in series) {

     call quality.cleanup as cleaning  { input:
        condition = condition
     }

     call quant.quantify as counting { input:
          transcriptome = transcriptomeIndex,
          condition = condition[0],
          samples = cleaning.out["trimmed"]
     }
    }

    output {
        Array[Pair[String, Array[File]]] out = counting.out
    }

}
