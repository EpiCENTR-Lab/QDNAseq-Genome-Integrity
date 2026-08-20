###############################################################################
# install_packages.R
#
# Installs all required packages for the QDNAseq genome integrity pipeline,
# pinned to the exact versions used to generate the published analysis
# (see SESSION_INFO.txt from that run).
#
#   R version:          4.5.3
#   Bioconductor release: 3.22
#
# Run once before using the pipeline:
#
#   Rscript install_packages.R
#
###############################################################################

message("Installing required packages for the QDNAseq pipeline...")

###############################################################################
# R version check
###############################################################################

required_r_version <- "4.5.3"
current_r_version <- paste(R.version$major, R.version$minor, sep = ".")

if (current_r_version != required_r_version) {
  message(
    "NOTE: this pipeline was built and validated on R ", required_r_version,
    ", but you are running R ", current_r_version, ". ",
    "Package versions below were pinned against R ", required_r_version,
    " — a different R version may not have binaries available for every ",
    "pinned version, or may resolve dependencies differently. Proceeding anyway."
  )
}

###############################################################################
# Install BiocManager and remotes (required for version pinning below)
###############################################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}

###############################################################################
# Pin the Bioconductor release
#
# Bioconductor packages are versioned together as part of a release; you
# cannot pin an individual Bioconductor package to an arbitrary version
# independent of the others. Pinning the release below reproduces the
# QDNAseq / Biobase / DNAcopy / Rsamtools / GenomeInfoDb versions used during
# pipeline development, since all of them were installed from this release.
###############################################################################

required_bioc_version <- "3.22"

BiocManager::install(version = required_bioc_version, ask = FALSE, update = FALSE)

###############################################################################
# CRAN packages, pinned to exact versions via remotes::install_version()
###############################################################################

cran_package_versions <- c(
  dplyr    = "1.2.1",
  tidyr    = "1.3.2",
  tibble   = "3.3.1",
  ggplot2  = "4.0.2",
  optparse = "1.7.5"
)

installed <- installed.packages()

for (pkg in names(cran_package_versions)) {
  wanted_version <- cran_package_versions[[pkg]]
  current_version <- if (pkg %in% rownames(installed)) {
    installed[pkg, "Version"]
  } else {
    NA_character_
  }

  if (is.na(current_version) || current_version != wanted_version) {
    message(
      "Installing ", pkg, " ", wanted_version,
      if (!is.na(current_version)) paste0(" (found ", current_version, ")") else ""
    )
    remotes::install_version(
      pkg,
      version = wanted_version,
      repos = "https://cloud.r-project.org",
      upgrade = "never"
    )
  } else {
    message(pkg, " ", wanted_version, " already installed.")
  }
}

###############################################################################
# Bioconductor packages
#
# Installed from the pinned Bioconductor release above, so these resolve to
# the exact versions used for the published analysis:
#   QDNAseq      1.46.0
#   Biobase      2.70.0
#   DNAcopy      1.84.0
#   Rsamtools    2.26.0
#   GenomeInfoDb (version not captured in SESSION_INFO.txt, since it loads
#                 as a dependency rather than being attached directly —
#                 pinned indirectly via the Bioconductor release)
###############################################################################

bioc_packages <- c(
  "QDNAseq",
  "Biobase",
  "DNAcopy",
  "Rsamtools",
  "GenomeInfoDb"
)

for (pkg in bioc_packages) {
  BiocManager::install(pkg, ask = FALSE, update = FALSE)
}

###############################################################################
# hg38 annotation package for QDNAseq
#
# QDNAseq.hg38 is not captured in sessionInfo()/SESSION_INFO.txt, since
# annotation packages are loaded via data() rather than library(). It is
# pinned here to the exact commit installed during pipeline development.
###############################################################################

qdnaseq_hg38_sha <- "cf7c07e39de0ac64a9c38cb030cba4626e2aae83"

installed_sha <- tryCatch(
  packageDescription("QDNAseq.hg38")$GithubSHA1,
  error = function(e) NULL,
  warning = function(w) NULL
)

if (is.null(installed_sha) || !identical(installed_sha, qdnaseq_hg38_sha)) {
  message(
    "Installing QDNAseq.hg38 @ ", substr(qdnaseq_hg38_sha, 1, 7),
    if (!is.null(installed_sha)) paste0(" (found ", substr(installed_sha, 1, 7), ")") else ""
  )
  remotes::install_github(paste0("asntech/QDNAseq.hg38@", qdnaseq_hg38_sha))
} else {
  message("QDNAseq.hg38 @ ", substr(qdnaseq_hg38_sha, 1, 7), " already installed.")
}

###############################################################################
# Verify installation
###############################################################################

message("\nVerifying installation...")

required <- c(names(cran_package_versions), bioc_packages, "QDNAseq.hg38")

missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]

if (length(missing) == 0) {
  message("All required packages installed successfully.")
} else {
  message("The following packages failed to install:")
  print(missing)
}

message("\nInstallation complete.")
