FROM continuumio/miniconda:4.7.12

WORKDIR /home/

COPY environment.yaml ./

RUN conda env create -f environment.yaml

RUN echo "source activate freebayes" > ~/.bashrc
ENV PATH /opt/conda/envs/freebayes/bin:$PATH

ENTRYPOINT ["snakemake"]