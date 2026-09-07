# base: miniforge ships conda + libmamba
FROM condaforge/miniforge3:latest

# git for cplfold clone
RUN apt-get update && apt-get install -y --no-install-recommends \
        git ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/hyb2
COPY . /opt/hyb2

# create env
RUN conda env create -f environment.yml && conda clean -afy

# put env on path
ENV PATH=/opt/conda/envs/hyb2/bin:$PATH \
    UNAFOLDDAT=/opt/conda/envs/hyb2/share/oligoarrayaux \
    HYB2_VARNA_JAR=/opt/hyb2/VARNA/build/jar/VARNAcmd.jar \
    HYB2_CPLFOLD_DIR=/opt/hyb2/CPLfold

# clone cplfold 
RUN git clone -b feature/standalone-pseudoknot-energy https://github.com/Vicky-0256/CPLfold.git /opt/hyb2/CPLfold

# check everything is in place
RUN sh /opt/hyb2/docker_scripts/preflight.sh

# ensure that the tmp_dir guard is in executable
RUN chmod +x /opt/hyb2/docker_scripts/entrypoint.sh
ENTRYPOINT ["/opt/hyb2/docker_scripts/entrypoint.sh"]
CMD ["hyb2", "--help"]
# trigger ci

