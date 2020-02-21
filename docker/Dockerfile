From debian:stretch

RUN apt-get update && \
  apt-get install -y build-essential \
    wget \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libncurses-dev \
    libbz2-dev \
    liblzma-dev \
    libxml2-dev \
    zip \
    python3-dev \
    python-dev \
    python3-pip \
    default-jre \
    r-base \
    pandoc \
    git

# Install all software under /opt/software:
RUN mkdir -p /opt/software 

# Install fastQC:
RUN cd /opt/software && \
  wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip && \
  unzip fastqc_v0.11.8.zip && \
  chmod +x /opt/software/FastQC/fastqc
ENV PATH="/opt/software/FastQC:${PATH}"

# Install HISAT2
RUN cd /opt/software && \
  wget https://cloud.biohpc.swmed.edu/index.php/s/hisat2-220-source/download && \
  unzip download && \
  cd hisat2-2.2.0/ && \
  make 
ENV PATH="/opt/software/hisat2-2.2.0:${PATH}"

# Install samtools:
RUN cd /opt/software && \
  wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
  tar -xjf samtools-1.9.tar.bz2 && \
  cd samtools-1.9 && \
  ./configure --prefix=/opt/software/samtools && \
  make && \
  make install
ENV PATH="/opt/software/samtools/bin:${PATH}"

# Get the Picard JAR:
RUN mkdir -p /opt/software/picard && \
  cd /opt/software/picard && \
  wget https://github.com/broadinstitute/picard/releases/download/2.18.26/picard.jar
ENV PICARD_JAR=/opt/software/picard/picard.jar

# subRead featureCounts
RUN cd /opt/software && \
  wget https://sourceforge.net/projects/subread/files/subread-1.6.3/subread-1.6.3-source.tar.gz && \
  tar -xzf subread-1.6.3-source.tar.gz && \
  cd subread-1.6.3-source/src && \
  make -f Makefile.Linux
ENV PATH="/opt/software/subread-1.6.3-source/bin:${PATH}"

# Install libStatGen and fastQValidator:
RUN cd /opt/software && \
  git clone https://github.com/statgen/libStatGen.git && \
  cd libStatGen && \
  make && \
  cd .. && \
  git clone https://github.com/statgen/fastQValidator.git && \
  cd fastQValidator && \
  make && \
  make install

# Install DESeq2 and dependencies:
RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('RCurl'); biocLite('XML'); biocLite('DESeq2')"
RUN R -e "install.packages('reshape2', repos = 'http://cran.us.r-project.org'); install.packages('tidyverse', repos = 'http://cran.us.r-project.org'); install.packages('ggdendro', repos = 'http://cran.us.r-project.org'); install.packages('gplots', repos = 'http://cran.us.r-project.org'); install.packages('ggplot2', repos = 'http://cran.us.r-project.org');"
ADD deseq2.R /opt/software/
ADD make_dge_plots.py /opt/software/
ADD make_figures.R /opt/software/
ADD draw_heatmap.R /opt/software/
ADD draw_pca.R /opt/software/
ADD contrast_independent_figures.R /opt/software/

# Install some Python3 libraries:
ADD requirements.txt /opt/software/
RUN pip3 install -r /opt/software/requirements.txt

# Now install RSeQC:
RUN pip3 install RSeQC

# Install multiqc:
RUN pip3 install multiqc

# Install Jinja2 for templating.  Likely installed with other python packages, but just to be sure:
RUN pip3 install Jinja2

# for multiqc to read commandline:
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

# Add other scripts under /usr/local/bin so they are in PATH:
# Script to concatenate the featureCounts output into a matrix:
ADD concatenate_featurecounts.py /usr/local/bin/
RUN chmod +x /usr/local/bin/concatenate_featurecounts.py

# The "alternate" script for inferring the strandedness:
ADD alternate_infer_experiment.py /usr/local/bin/
RUN chmod +x /usr/local/bin/alternate_infer_experiment.py

# The scripts for performing pre-workflow checks:
RUN mkdir /opt/software/precheck
ADD perform_precheck.py /opt/software/precheck/
ADD check_fastq.py /opt/software/precheck/

# The script for generating the markdown report:
ADD generate_report.py /usr/local/bin/
RUN chmod +x /usr/local/bin/generate_report.py

# The script for reformatting the annotations:
ADD reformat_annotations.py /opt/software/

# Add the report template and styling:
RUN mkdir /opt/report
ADD report.md /opt/report/
ADD report.css /opt/report/

ENTRYPOINT ["/bin/bash"]

