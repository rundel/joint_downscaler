packages = installed.packages()[,1]

needed_packages = c("plyr","lubridate","Rcpp","RcppArmadillo","raster","fields","coda")

for(p in needed_packages)
{
    if (!(p %in% packages))
        install.packages(p)
}

#update.packages(needed_packages)
#installed.packages()[installed.packages()[,1] %in% needed_packages,3]
