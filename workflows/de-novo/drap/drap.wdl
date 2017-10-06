workflow Drap {
# https://peerj.com/articles/2988/
  File reads_1
  File reads_2
  File proteins #homology proteins to check assembly


}



task runDrap {
#runDrap performs the assembly including compaction and correction.
#It produces a contig set but also a HTML log report presenting different assembly metrics

}

task runMeta {
#runMeta merges and compacts different contigs sets
#and should be used for very large datasets for which memory
#or CPU requirements do not enable a unique global assembly or for highly complex datasets

}

task runAssessment {
#runAssessment compares different contig sets and gathers the results in a global report.

}