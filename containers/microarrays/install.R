url <- "http://bioconductor.org/packages/3.6/bioc"

if ("BiocInstaller" %in% rownames(installed.packages()))
remove.packages("BiocInstaller")
install.packages("BiocInstaller", repos=url)

to_install <- c("Matrix", "KernSmooth", "mgcv", "devtools", "biomaRt", "limma",
                "affy", "lumi", "methylumi", "minfi", "stringr", "GEOquery", "GEOmetadb",
                "ExiMiR", "AgiMicroRna", "doParallel","foreach", "RJSONIO",
                "IlluminaHumanMethylation27k.db", "FDb.InfiniumMethylation.hg19"
)

for (pack in to_install)
    if (!suppressWarnings(require(pack, character.only=TRUE)))
        BiocInstaller::biocLite(pack)
suppressWarnings(BiocInstaller::biocValid(fix=TRUE, ask=FALSE))
suppressWarnings(install.packages("optparse"))