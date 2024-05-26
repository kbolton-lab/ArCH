#!/usr/bin/env Rscript

packages = commandArgs(trailingOnly=TRUE)

for (l in packages) {

    install.packages(l, dependencies=TRUE, repos='https://cran.rstudio.com/');

    if ( ! library(l, character.only=TRUE, logical.return=TRUE) ) {
        quit(status=1, save='no')
    }
}