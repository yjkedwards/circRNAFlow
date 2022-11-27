if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos = "http://cran.us.r-project.org",version="3.15")
#install it here
BiocManager::install("clusterProfiler",verion="4.4.4",force=TRUE)

