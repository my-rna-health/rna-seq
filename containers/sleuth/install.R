url <- "http://bioconductor.org/packages/3.5/bioc"

if ("BiocInstaller" %in% rownames(installed.packages()))
remove.packages("BiocInstaller")
install.packages("BiocInstaller", repos=url)

to_install <- c("Matrix", "KernSmooth", "mgcv", "devtools", "pachterlab/sleuth", "COMBINE-lab/wasabi", "biomaRt")
for (pack in to_install)
    if (!suppressWarnings(require(pack, character.only=TRUE)))
        BiocInstaller::biocLite(pack)
suppressWarnings(BiocInstaller::biocValid(fix=TRUE, ask=FALSE))
suppressWarnings(install.packages("optparse"))
