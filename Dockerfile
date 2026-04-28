FROM continuumio/miniconda3:latest

LABEL maintainer="Brown Beckley <brownbeckley94@gmail.com>"
LABEL description="AcinetoScope - Comprehensive Acinetobacter baumannii Genomic Typing Platform"

# Install system dependencies (procps for resource monitoring, jq for JSON parsing)
RUN apt-get update && apt-get install -y \
    procps \
    jq \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /opt/acinetoscope

# Copy entire project
COPY . /opt/acinetoscope/

# Create the Conda environment from environment.yml
RUN conda env create -f environment.yml && \
    conda clean -afy

# Make the environment the default for RUN commands
SHELL ["conda", "run", "-n", "acinetoscope", "/bin/bash", "-c"]

# Run abricate database setup (one-time)
RUN abricate --setupdb

# Set entrypoint to acinetoscope command
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "acinetoscope", "acinetoscope"]
CMD ["-h"]
