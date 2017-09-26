#https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html

suppressMessages({
    library("sleuth")
})
series <- "results"
conditions <- "hiseq_info.txt"
sample_id <- dir(file.path(series))
kal_dirs <- file.path(sample_id, "kallisto")

s2c <- read.table(file.path("data", "metadata", conditions), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample = run_accession, condition)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
print(s2c)

# the sleuth object must first be initialized with
so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)

# then the full model is fit with
so <- sleuth_fit(so, ~condition, 'full')

# the “reduced” model is fit with
so <- sleuth_fit(so, ~1, 'reduced')

# and the test is performed with
so <- sleuth_lrt(so, 'reduced', 'full')
models(so)

# the results of the test can be examined with
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
head(sleuth_significant, 20)

# the table shown above displays the top 20 significant genes with a (Benjamini-Hochberg multiple testing corrected) q-value <= 0.05.
plot_bootstrap(so, "ENST00000263734", units = "est_counts", color_by = "condition")