# SeqGuide: one image with the app and every tool it needs.
# Builds natively on Intel/AMD (Windows PCs) and on Apple silicon Macs.
FROM mambaorg/micromamba:latest

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes

# Medaka (Nanopore polishing) has heavy dependencies, so it lives in its own environment.
RUN micromamba create -y -n medaka -c conda-forge -c bioconda medaka=2.2.2 && micromamba clean --all --yes

USER root
RUN apt-get update && apt-get install -y --no-install-recommends tzdata && rm -rf /var/lib/apt/lists/*
RUN printf '#!/bin/bash\nexec micromamba run -n medaka medaka_consensus "$@"\n' > /usr/local/bin/medaka_consensus && \
    printf '#!/bin/bash\nexec micromamba run -n medaka medaka "$@"\n' > /usr/local/bin/medaka && \
    chmod +x /usr/local/bin/medaka_consensus /usr/local/bin/medaka && \
    mkdir -p /data && chown $MAMBA_USER:$MAMBA_USER /data
USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER seqguide /opt/seqguide/seqguide
COPY --chown=$MAMBA_USER:$MAMBA_USER app /opt/seqguide/app
COPY --chown=$MAMBA_USER:$MAMBA_USER .streamlit /opt/seqguide/.streamlit

ENV SEQGUIDE_DATA=/data
WORKDIR /opt/seqguide
EXPOSE 8501
CMD ["streamlit", "run", "app/streamlit_app.py", "--server.address=0.0.0.0", "--server.port=8501"]
