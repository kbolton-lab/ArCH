FROM bioconductor/bioconductor_docker:RELEASE_3_14

ADD install_packages_or_die.R /
ADD install_bioconductor_or_die.R /
ADD install_devtools_or_die.R /

RUN Rscript --no-save install_packages_or_die.R tidyverse sqldf jsonlite pryr
RUN Rscript --no-save install_packages_or_die.R optparse vcfR dplyr stringr
Run Rscript --no-save install_bioconductor_or_die.R GenomeInfoDb
RUN Rscript --no-save install_bioconductor_or_die.R BSgenome BSgenome.Hsapiens.UCSC.hg38
RUN Rscript --no-save install_packages_or_die.R tidyselect
RUN Rscript --no-save install_bioconductor_or_die.R R453Plus1Toolbox
RUN Rscript --no-save install_bioconductor_or_die.R Biostrings
RUN Rscript --no-save install_packages_or_die.R

RUN mkdir /opt/bin
RUN apt-get update && apt-get -y install less vim bcftools

RUN apt-get clean && \
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
apt-get autoclean && \
apt-get autoremove -y && \
rm -rf /var/lib/{apt,dpkg,cache,log}/

RUN mkdir /opt/bin/annotate/
RUN rm -rf /opt/bin/annotate/ArCHAnnotationScript.R && rm -rf /opt/bin/annotate/supportFunctions.R
COPY ./annotate/ArCHAnnotationScript.R /opt/bin/annotate/ArCHAnnotationScript.R
COPY ./annotate/supportFunctions.R /opt/bin/annotate/supportFunctions.R

RUN mkdir /opt/bin/combine/
RUN rm -rf /opt/bin/combine/ArCHCombine.R && rm -rf /opt/bin/combine/supportFunctions.R
COPY ./combine/ArCHCombine.R /opt/bin/combine/ArCHCombine.R
COPY ./combine/supportFunctions.R /opt/bin/combine/supportFunctions.R

RUN mkdir /opt/bin/post_pipeline/
RUN  rm -rf /opt/bin/post_pipeline/ArCHPostPipeline.R && rm -rf /opt/bin/post_pipeline/supportFunctions.R
COPY ./post_pipeline/ArCHPostPipeline.R /opt/bin/post_pipeline/ArCHPostPipeline.R
COPY ./post_pipeline/supportFunctions.R /opt/bin/post_pipeline/supportFunctions.R

CMD ["/bin/bash"]
