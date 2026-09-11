FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/13502_wf_init_archrproject_workflow:0.27.2-cc2abe-wip-9c8f70

RUN mkdir /opt/latch

# Install specific version of numpy
RUN pip install numpy==1.25.2

# Install pip dependencies from `requirements.txt`, pychromvar with chunks
COPY requirements.txt /opt/latch/requirements.txt
RUN pip install --requirement /opt/latch/requirements.txt

RUN pip3 uninstall -y aiobotocore botocore awscli s3transfer
RUN pip3 install awscli

RUN R -e "remotes::install_github('jpmcga/ArchR', ref = '619f75d')"
RUN R -e "BiocManager::install(c('BSgenome.Mmusculus.UCSC.mm39', 'TxDb.Mmusculus.UCSC.mm39.knownGene', 'org.Mm.eg.db'), ask = FALSE, update = FALSE)"

# Pin compatible plotting packages after other R installs, which can upgrade them.
RUN R -e 'remotes::install_version("ggplot2", version = "3.4.1", repos = "https://cran.r-project.org", upgrade = "never"); remotes::install_version("ggrepel", version = "0.9.6", repos = "https://cran.r-project.org", upgrade = "never")'
# Fail the image build if installation failed or these namespaces cannot coexist.
RUN R -e 'stopifnot(packageVersion("ggplot2") == "3.4.1", packageVersion("ggrepel") == "0.9.6"); library(ggplot2); library(ggrepel); library(ArchR)'

# Copy output files for Shiny app
COPY getDeviation_ArchR.R /root/getDeviation_ArchR.R

COPY archrproject.Rproj /root/archrproject.Rproj
COPY .renvignore /root/.renvignore

COPY custom_ArchR_genomes_and_annotations /root/custom_ArchR_genomes_and_annotations

# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
RUN python3 -m pip install latch==2.66.3
RUN pip install https://github.com/atlasxomics/atx-common/archive/refs/tags/v0.1.0.tar.gz
RUN python3 -m pip install setuptools==61.0.0

RUN rm -r /root/wf
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root
