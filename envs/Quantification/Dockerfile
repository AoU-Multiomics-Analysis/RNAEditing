FROM mambaorg/micromamba:1.5.3

ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV PERLLIB=/opt/scripts   

RUN micromamba install -y -n base -c conda-forge -c bioconda \
    perl \
    perl-io-compress \
    python=3.11.1 \
    bioconda::samtools=1.23 \
    bioconda::htslib=1.23 \
    conda-forge::numpy \
    conda-forge::scipy \
    conda-forge::pandas \
    conda-forge::r-base \
    conda-forge::r-optparse \
    conda-forge::r-tidyverse \
    conda-forge::r-data.table \
    bioconda::bioconductor-pcatools \
    && micromamba clean --all --yes

WORKDIR /opt/scripts

COPY ../quantification/1_pileup/scripts/query_editing_levels.pl /opt/scripts/
COPY ../quantification/1_pileup/scripts/parse_pileup_query.pl /opt/scripts/
COPY ../quantification/2_combine_matrices/scripts/combine_sample_matrices.pl /opt/scripts/
COPY ../quantification/3_transform/scripts/transform.py /opt/scripts/
COPY ../quantification/3_transform/scripts/calculate_pcs.R /opt/scripts/


COPY ../quantification/1_pileup/references/augmented_editing_sites_REDIportal_plus_original_annotated.bed.gz /opt/scripts/

ENV PATH="/opt/scripts:${PATH}"

CMD ["bash"]
