#!/bin/bash

#
# Copyright 2021 Simone Maestri. All rights reserved.
# Simone Maestri <simone.maestri@univr.it>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

#number of threads
THREADS=8
#ploidy of the sample
PLOIDY=2
########################################################################################################
PIPELINE_DIR="/path/to/Illumina_variant_calling"
MINICONDA_DIR="/path/to/miniconda3"
########### End of user editable region #################################################################
FASTQC=$MINICONDA_DIR"/envs/Illumina_variant_calling_env/bin/fastqc"
FASTP=$MINICONDA_DIR"/envs/Illumina_variant_calling_env/bin/fastp"
BWA=$MINICONDA_DIR"/envs/Illumina_variant_calling_env/bin/bwa"
SAMTOOLS=$MINICONDA_DIR"/envs/Illumina_variant_calling_env/bin/samtools"
JAVA=$MINICONDA_DIR"/envs/Illumina_variant_calling_env/bin/java"
PICARD=$MINICONDA_DIR"/envs/Illumina_variant_calling_env/bin/picard"
FGBIO=$MINICONDA_DIR"/envs/Illumina_variant_calling_env/bin/fgbio"
GATK=$PIPELINE_DIR_DIR"/gatk-4.2.0.0/gatk"
QUALIMAP=$MINICONDA_DIR"/envs/Illumina_variant_calling_env/bin/qualimap"
