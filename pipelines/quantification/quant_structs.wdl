version development

struct QuantifiedRun {
    String run
    File folder
    File quant
    File lib
    #Map[String, String] metadata
    Array[Pair[String, String]] metainfo
}

struct QuantifiedGSM {
    Array[QuantifiedRun] runs
    File metadata
}