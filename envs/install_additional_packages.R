# Run with `Rscript envs/install_additional_packages.R`
# Install spacexr directly from Github
# Skip dependencies packages updates
devtools::install_github("dmcable/spacexr",
                         build_vignettes = FALSE,
                         upgrade = "never")