##FROM mambaorg/micromamba:1.5.3
#FROM ghcr.io/prefix-dev/pixi:latest

### Set up environment
##ENV MAMBA_DOCKERFILE_ACTIVATE=1
##ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

### Create a new environment and install packages
#RUN pixi global install -c conda-forge -c bioconda python=3.11.1  google-cloud-sdk bioconda::samtools=1.23 perl bioconda::htslib=1.23
#RUN pixi global install --environment google-cloud-sdk python=3.11 
#RUN pixi global install --environment python \
    #conda-forge::numpy \
    #conda-forge::scipy \
    #conda-forge::pandas 
FROM ghcr.io/prefix-dev/pixi:latest

RUN pixi global install --environment python \
    -c conda-forge -c bioconda \
    python=3.11.1 \
    google-cloud-sdk \
    bioconda::samtools=1.23 \
    perl \
    bioconda::htslib=1.23

RUN pixi global add --environment python \
    conda-forge::numpy \
    conda-forge::scipy \
    conda-forge::pandas


COPY /quantification/1_pileup/scripts/query_editing_levels.pl /opt/scripts/
COPY quantification/1_pileup/scripts/parse_pileup_query.pl /opt/scripts/
COPY /quantification/2_combine_matrices/scripts/combine_sample_matrices.pl /opt/scripts/
COPY quantification/3_transform/scripts/transform.py /opt/scripts/

COPY /quantification/1_pileup/references/All.AG.stranded.annovar.Hg38_multianno.AnnoAlu.AnnoRep.NR.bed.gz /opt/scripts/

