if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos = "http://cran.us.r-project.org")

BiocManager::install("biomaRt")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")


