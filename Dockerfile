FROM mambaorg/micromamba:1.5.3

ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

RUN micromamba install -y -n base -c conda-forge -c bioconda \
    perl=5.34 \
    python=3.11.1 \
    bioconda::samtools=1.23 \
    bioconda::htslib=1.23 \
    conda-forge::numpy \
    conda-forge::scipy \
    conda-forge::pandas \
    conda-forge::gzip \
    && micromamba clean --all --yes

RUN mkdir -p /code
WORKDIR /code

COPY quantification/1_pileup/scripts/query_editing_levels.pl /code/
COPY quantification/1_pileup/scripts/parse_pileup_query.pl /code/
COPY quantification/2_combine_matrices/scripts/combine_sample_matrices.pl /code/
COPY quantification/3_transform/scripts/transform.py /code/

COPY quantification/1_pileup/references/All.AG.stranded.annovar.Hg38_multianno.AnnoAlu.AnnoRep.NR.bed.gz /code/

RUN chmod +x /code/*.pl /code/*.py

ENV PATH="/code:${PATH}"

CMD ["bash"]