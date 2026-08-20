# base: miniforge ships conda + libmamba
FROM condaforge/miniforge3:latest

# builds deps for hotknots + git for cplfold clone
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential git ca-certificates libx11-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/hyb2
COPY . /opt/hyb2

# create env
RUN conda env create -f environment.yml && conda clean -afy

# put env on path so hyb2-py resolves
ENV PATH=/opt/conda/envs/hyb2/bin:$PATH \
    UNAFOLDDAT=/opt/conda/envs/hyb2/share/oligoarrayaux \
    HYB2_VARNA_JAR=/opt/hyb2/VARNA/build/jar/VARNAcmd.jar \
    HYB2_CPLFOLD_DIR=/opt/hyb2/CPLfold

# clone cplfold + build hotknots
RUN git clone https://github.com/Vicky-0256/CPLfold.git /opt/hyb2/CPLfold \
    && cd /opt/hyb2/CPLfold \
    && git checkout 6f49167b25cf6312c45410a4a3b7f22ae1ca54b9 \
    && cd Utils/HotKnots_v2.0 \
    && find . -name '*.o' -delete && find . -name '*.a' -delete \
    && cd graphics && make graphics.o PlotRna.o \
    && cd .. && make

# check everything is in place
RUN sh /opt/hyb2/docker_scripts/preflight.sh

# ensure that the tmp_dir guard is in executable
RUN chmod +x /opt/hyb2/docker_scripts/entrypoint.sh
ENTRYPOINT ["/opt/hyb2/docker_scripts/entrypoint.sh"]
CMD ["hyb2-py", "--help"]