FROM mambaorg/micromamba:jammy
COPY --chown=$MAMBA_USER:$MAMBA_USER peakpicking.R /osfd/peakpicking.R
COPY --chown=$MAMBA_USER:$MAMBA_USER peakpicking2.cpp /osfd/peakpicking2.cpp
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml /tmp/env.yaml

USER root
RUN apt-get update && apt-get install -y build-essential
RUN mkdir /data/

USER $MAMBA_USER
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes
