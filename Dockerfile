FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:fe0b-main
RUN apt-get update -y
RUN apt-get install -y gdebi-core 
RUN apt install -y aptitude
RUN aptitude install -y libjpeg-dev

# Install R
RUN apt-get update -y && \
    apt-get install -y \
        r-base \
        r-base-dev \
        apt-transport-https \
        build-essential \
        gfortran \
        libhdf5-dev \
        libatlas-base-dev \
        libbz2-dev \        
        libcurl4-openssl-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libgdal-dev \
        libgit2-dev \
        libgsl-dev \
        libicu-dev \
        liblzma-dev \
        libmagick++-dev \
        libpango-1.0-0 \
        libpangocairo-1.0-0 \
        libpcre3-dev \
        libssl-dev \
        libtcl8.6 \
        libtiff5 \
        libtk8.6 \
        libxml2-dev \
        libxt-dev \
        libx11-dev \
        libtiff-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        locales \
        make \
        pandoc \
        r-cran-rjava \
        tzdata \
        vim \
        wget \
        zlib1g-dev

# Upgrade R to version 4.3.0
RUN wget https://cran.r-project.org/src/base/R-4/R-4.3.0.tar.gz
RUN tar zxvf R-4.3.0.tar.gz
RUN cd R-4.3.0 && ./configure --enable-R-shlib
RUN cd R-4.3.0 && make && make install

# Have to install devtools, cairo like this; see https://stackoverflow.com/questions/20923209
RUN apt-get install -y r-cran-devtools libcairo2-dev

# Install java
RUN apt install -y default-jdk
RUN R CMD javareconf

# Fix systemd conflict with timedatectl
RUN echo "TZ=$( cat /etc/timezone )" >> /etc/R/Renviron.site

# Installation of R packages with renv
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/renv_1.0.5.tar.gz', repos = NULL, type = 'source')"
COPY renv.lock /root/renv.lock
RUN mkdir /root/renv
COPY renv/activate.R /root/renv/activate.R
RUN R -e "renv::restore()"

# Need to specify the order, numpy first
RUN python3 -m pip install numpy==1.26.2
RUN python3 -m pip install MACS2==2.2.6

# Copy output files for Shiny app
COPY getDeviation_ArchR.R /root/getDeviation_ArchR.R
COPY makeShinyFiles.R /root/makeShinyFiles.R

COPY server.R /root/server.R
COPY ui.R /root/ui.R

COPY server_2.R /root/server_2.R
COPY ui_2.R /root/ui_2.R

COPY server_3.R /root/server_3.R
COPY ui_3.R /root/ui_3.R

COPY server_4.R /root/server_4.R
COPY ui_4.R /root/ui_4.R

COPY www /root/www

# Latch SDK
# DO NOT REMOVE
RUN pip install latch==2.40.2

RUN mkdir /opt/latch

# Latch workflow registration metadata
# DO NOT CHANGE
ARG tag
# DO NOT CHANGE
ENV FLYTE_INTERNAL_IMAGE $tag

COPY wf /root/wf
WORKDIR /root