# Use the Bioconductor Docker image as the base image
FROM bioconductor/bioconductor_docker:RELEASE_3_19

# Make a directory in the container
RUN mkdir /home/GENOXA2

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*
    
# Install LaTeX
RUN apt-get update && apt-get install -y \
    texlive \
    texlive-latex-extra \
    texlive-fonts-recommended \
    && rm -rf /var/lib/apt/lists/*

# Install required R packages from CRAN
RUN R -e "install.packages(c('parallel','parallelly', 'doParallel','stats', 'stats4', 'tidyr','ggrepel','viridis','viridisLite', 'grid','rrcovHD','naivebayes','elasticnet','msaenet', 'httpuv', 'jsonlite', 'shiny', 'shinyjs', 'shinydashboard', 'DT', 'caret', 'grid', 'naivebayes', 'reshape2', 'pROC', 'shinyalert', 'tidyverse', 'ggplot2', 'dplyr', 'pheatmap', 'shinycssloaders', 'ggupset', 'cowplot', 'gridExtra', 'msaenet', 'keras', 'xgboost', 'e1071', 'mboost', 'mda', 'plyr', 'rotationForest', 'C50', 'earth', 'ipred', 'party', 'sdwd', 'kernlab', 'glmnet', 'randomForest', 'kerndwd','MASS','ordinalForest','partykit','ranger','rrcovHD','ada','adabag','arm','bartMachine','binda','bnclassify', 'bst','caTools','class','frbs','gam','hda','HDclassif','HiDimDA','klaR','LiblineaR','mgcv', 'nnet','ordinalNet','pamr','proxy','randomGLM','rFerns','robustDA','rocc','rrcov','rrcovHD', 'RWeka','sda','snn','sparsediscrim','sparseLDA','stepPlr','supervisedPRIM','VGAM','wsrf'))"

# Install Bioconductor packages
RUN R -e "BiocManager::install(c('S4Vectors','IRanges','Biobase', 'BiocGenerics', 'limma', 'pcaMethods', 'edgeR', 'clusterProfiler', 'pathview', 'AnnotationDbi', 'enrichplot', 'org.Hs.eg.db', 'org.Rn.eg.db', 'org.Sc.sgd.db', 'org.Dm.eg.db', 'org.At.tair.db', 'org.Dr.eg.db', 'org.Bt.eg.db', 'org.Ce.eg.db', 'org.Gg.eg.db', 'org.Ss.eg.db', 'org.Mm.eg.db', 'org.EcK12.eg.db', 'org.Cf.eg.db', 'org.Xl.eg.db', 'org.Pt.eg.db', 'org.Ag.eg.db', 'org.EcSakai.eg.db'))"

# Copy your Shiny app code to the image
COPY . /home/GENOXA2

# Set working directory
WORKDIR /home/GENOXA2

# Set an environment variable to prevent xdg-open errors
ENV BROWSER=none

# Expose the port for Shiny
EXPOSE 80

# Run the Shiny app
CMD ["R", "-e", "shiny::runApp('/home/GENOXA2', host='0.0.0.0', port=80)"]