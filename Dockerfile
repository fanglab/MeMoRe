# Base image https://hub.docker.com/u/rocker/
FROM rocker/shiny:4.0.4

# Install debian packages
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    libxml2-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libpq-dev \
    libssh2-1-dev \
    unixodbc-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    git

# Update system libraries
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean

# Retrieve repository
RUN git clone https://github.com/fanglab/MeMoRe.git

# Define WORKDIR and flag as Docker container
WORKDIR /MeMoRe
RUN touch iam_a_container

# Install renv & restore packages
RUN Rscript -e 'install.packages("renv")'
RUN Rscript -e 'renv::consent(provided = TRUE)'
RUN Rscript -e 'renv::restore()'

# Expose port
EXPOSE 3838

# Run shiny app on container start
CMD ["R", "-e", "shiny::runApp('./app', host = '0.0.0.0', port = 3838)"]
