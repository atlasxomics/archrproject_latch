FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/13502_wf_init_archrproject_workflow:0.17.1-d3928b
RUN apt-get update -y

COPY getDeviation_ArchR.R /root/getDeviation_ArchR.R
COPY makeShinyFiles.R /root/makeShinyFiles.R
COPY www /root/www

COPY uiserver50by50 /root/uiserver50by50
COPY uiserver96by96 /root/uiserver96by96
# numpy needed to be install before macs2 v-2.2.6
COPY requirements.txt /root/requirements.txt
RUN python3 -m pip install -r requirements.txt

# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
#RUN python3 -m pip install --upgrade latch
RUN python3 -m pip install latch==2.40.2

COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root
