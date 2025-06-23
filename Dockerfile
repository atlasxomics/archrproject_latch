FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/13502_wf_init_snap_workflow:0.17.25.6-5cf1a1-wip-66bb71

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
RUN python3 -m pip install --upgrade latch

RUN rm -r /root/wf
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root
