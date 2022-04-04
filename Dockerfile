FROM continuumio/miniconda3

########### set variables
ENV DEBIAN_FRONTEND noninteractive

########## generate working directories
RUN mkdir /home/tools

######### dependencies
RUN apt-get update -qq \
    && apt-get install -y \
    build-essential \
    wget \
    unzip \
    bzip2 \
    git \
    libidn11* \
    nano \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

############################################################ install Illumina_variant_calling
WORKDIR /home/tools/

RUN git clone https://github.com/MaestSi/Illumina_variant_calling.git
WORKDIR /home/tools/Illumina_variant_calling
RUN chmod 755 *

RUN sed -i 's/PIPELINE_DIR=.*/PIPELINE_DIR=\"\/home\/tools\/Illumina_variant_calling\/\"/' config_Variant_calling_pipeline.sh
RUN sed -i 's/MINICONDA_DIR=.*/MINICONDA_DIR=\"\/opt\/conda\/\"/' config_Variant_calling_pipeline.sh

RUN conda config --add channels r && \
conda config --add channels anaconda && \
conda config --add channels conda-forge && \
conda config --add channels bioconda

RUN conda create -n Illumina_variant_calling_env fastqc fastp fgbio picard bwa samtools=1.15 java-jdk qualimap
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.2.0.0/gatk-4.2.0.0.zip && unzip gatk-4.2.0.0.zip

WORKDIR /home/
