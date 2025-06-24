FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/13502_wf_init_archrproject_workflow:0.27.2-cc2abe-wip-9c8f70

RUN pip install latch==2.52.2
RUN mkdir /opt/latch

# Install specific version of numpy
RUN pip install numpy==1.25.2

# Install pip dependencies from `requirements.txt`, pychromvar with chunks
COPY requirements.txt /opt/latch/requirements.txt
RUN pip install --requirement /opt/latch/requirements.txt

RUN pip3 uninstall -y aiobotocore botocore awscli s3transfer
RUN pip3 install awscli

# Copy output files for Shiny app
COPY getDeviation_ArchR.R /root/getDeviation_ArchR.R
COPY makeShinyFiles.R /root/makeShinyFiles.R
COPY www /root/www

COPY archrproject.Rproj /root/archrproject.Rproj
COPY .renvignore /root/.renvignore

COPY uiserver50by50 /root/uiserver50by50
COPY uiserver96by96 /root/uiserver96by96

COPY custom_ArchR_genomes_and_annotations /root/custom_ArchR_genomes_and_annotations

# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
# RUN python3 -m pip install --upgrade latch

RUN rm -r /root/wf
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root
