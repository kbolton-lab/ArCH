#!/usr/bin/env Rscript

packages = commandArgs(trailingOnly=TRUE)

for (l in packages) {

    devtools:::install_github(l, dependencies=TRUE);

    if ( ! library(basename(l), character.only=TRUE, logical.return=TRUE) ) {
        quit(status=1, save='no')
    }
}