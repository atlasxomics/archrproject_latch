FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:dd8f-main


# Install R
RUN apt-get update -y && \
    apt-get install -y software-properties-common && \
    add-apt-repository "deb http://cloud.r-project.org/bin/linux/debian buster-cran40/" && \
    apt-get install -y r-base r-base-dev libxml2-dev libcurl4-openssl-dev libssl-dev wget

# Install packages
RUN R -e "install.packages('devtools')"
RUN R -e "install.packages('BiocManager')"
RUN R -e "devtools::install_github('GreenleafLab/ArchR', ref='master', repos = BiocManager::repositories())"
RUN R -e "library(ArchR); ArchR::installExtraPackages()"

RUN R -e "BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")"
RUN R -e "BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")"


# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
RUN python3 -m pip install --upgrade latch
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root

# RUN apt install -y dirmngr apt-transport-https ca-certificates software-properties-common gnupg2
# RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7'
# RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/debian buster-cran40/'
# RUN apt update
# RUN apt install -y r-base build-essential
# RUN apt-get install libcurl4-openssl-dev
