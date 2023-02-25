workflow quantification {
    input {
        Map[String, Directory] salmon_indexes
        Map[String, File] transcripts2genes
    }

    output {
        Array[MappedRun] mapped_runs = mapped_run
    }
}