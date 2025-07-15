FROM rocker/rstudio:4.4
SHELL ["/bin/bash", "-c"] 
RUN echo "" >/etc/apt/sources.list.d/ubuntu.sources && \
    echo "Types: deb" >>/etc/apt/sources.list.d/ubuntu.sources && \
    echo "URIs: https://mirrors.tuna.tsinghua.edu.cn/ubuntu/" >>/etc/apt/sources.list.d/ubuntu.sources && \
    echo "Suites: noble noble-updates noble-backports" >>/etc/apt/sources.list.d/ubuntu.sources && \
    echo "Components: main universe restricted multiverse" >>/etc/apt/sources.list.d/ubuntu.sources && \
    echo "Signed-By: /usr/share/keyrings/ubuntu-archive-keyring.gpg" >>/etc/apt/sources.list.d/ubuntu.sources && \
    echo "" >>/etc/apt/sources.list.d/ubuntu.sources && \
    echo "Types: deb" >>/etc/apt/sources.list.d/ubuntu.sources && \
    echo "URIs: https://mirrors.tuna.tsinghua.edu.cn/ubuntu/" >>/etc/apt/sources.list.d/ubuntu.sources && \
    echo "Suites: noble-security" >>/etc/apt/sources.list.d/ubuntu.sources && \
    echo "Components: main universe restricted multiverse" >>/etc/apt/sources.list.d/ubuntu.sources && \
    echo "Signed-By: /usr/share/keyrings/ubuntu-archive-keyring.gpg" >>/etc/apt/sources.list.d/ubuntu.sources && \
    #echo "" > /etc/apt/sources.list && \
    #echo "deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ jammy main restricted universe multiverse" >> /etc/apt/sources.list && \
    #echo "deb-src https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ jammy main restricted universe multiverse" >> /etc/apt/sources.list && \
    #echo "deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ jammy-updates main restricted universe multiverse" >> /etc/apt/sources.list && \
    #echo "deb-src https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ jammy-updates main restricted universe multiverse" >> /etc/apt/sources.list && \
    #echo "deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ jammy-backports main restricted universe multiverse" >> /etc/apt/sources.list && \
    #echo "deb-src https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ jammy-backports main restricted universe multiverse" >> /etc/apt/sources.list && \
    #echo "deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ jammy-security main restricted universe multiverse" >> /etc/apt/sources.list && \
    #echo "deb-src https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ jammy-security main restricted universe multiverse" >> /etc/apt/sources.list && \


    apt-get clean && \
    apt-get update && \
    apt-get install -y aptitude && \
    apt-get install -y pip && \
    apt-get install -y openjdk-8-jdk && \
    apt-get install -y zlib1g-dev && \
    apt-get install -y libhdf5-dev && \
    apt-get install -y libtirpc-dev && \
    apt-get install -y libbz2-dev && \
    apt-get install -y liblzma-dev && \
    #apt-get install -y libtinfo-dev=6.3-2ubuntu0.1 && \
    apt remove libncurses-dev libncurses5-dev libncursesw5-dev && \
    apt-get install -y libncurses5-dev libncursesw5-dev && \
    apt-get install -y libgsl-dev && \
    apt-get install -y cmake && \
    apt-get install -y libboost-iostreams-dev && \
#    apt-get install -y gcc-12 && \
#    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-12 100 && \
    apt-get install -y gcc-multilib && \
	ln /lib/x86_64-linux-gnu/libboost_iostreams.so.1.83.0 /lib/x86_64-linux-gnu/libboost_iostreams.so.1.71.0 && \
    apt-get install -y bedtools && \
    apt-get install -y libglpk40 libglpk-dev && \
    cp /bin/python3 /bin/python && \
    #pip install python3-pandas python3-numpy python3-scipy python3-scikit-learn && \
    #apt-get install -y python3-pandas && \
    #apt-get install -y python3-numpy && \
    #apt-get install -y python3-scipy && \
    #apt-get install -y python3-sklearn && \
    mkdir /software && \

    # install conda
    mkdir /software/miniconda && \
    cd /software/miniconda && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -u -p /software/miniconda3 && \
    rm /software/miniconda -rf
ENV PATH=$PATH:/software/miniconda3/bin
RUN cd /software && \
    wget https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 && \
    bzip2 -d htslib-1.21.tar.bz2 && \
    tar -xvf htslib-1.21.tar && \
    cd htslib-1.21 && \
    ./configure --prefix=/software/htslib/ && \
    make -j 8 && \
    make install   && \  
    rm ../htslib-1.21.tar ../htslib-1.21 -rf && \
    # install featureCoutns
    cd /software && \
    wget https://sourceforge.net/projects/subread/files/subread-2.0.6/subread-2.0.6-source.tar.gz && \
    tar -zxvf subread-2.0.6-source.tar.gz && \
    cd subread-2.0.6-source/src && \
    make -f Makefile.Linux -j 8 && \
    # install samtools
    cd /software && \
    wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
    bzip2 -d samtools-1.21.tar.bz2 && \
    tar -xvf samtools-1.21.tar && \
    cd samtools-1.21 && \
    ./configure --prefix=/software/samtools && \
    make -j 8 && \
    make install && \
    rm /software/samtools-1.21.tar && \
    conda create -n MAJIQ python=3.11 -y && \
    conda create -n SCSES python=3.11 -y
RUN cd /software && \
    wget https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 && \
    bzip2 -d htslib-1.21.tar.bz2 && \
    tar -xvf htslib-1.21.tar && \
    cd htslib-1.21 && \
    ./configure --prefix=/software/htslib/ && \
    make -j 8 && \
    make install   && \  
    rm ../htslib-1.21.tar ../htslib-1.21 -rf && \
    # install featureCoutns
    cd /software && \
    wget https://sourceforge.net/projects/subread/files/subread-2.0.6/subread-2.0.6-source.tar.gz && \
    tar -zxvf subread-2.0.6-source.tar.gz && \
    cd subread-2.0.6-source/src && \
    make -f Makefile.Linux -j 8 && \
    # install samtools
    cd /software && \
    wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
    bzip2 -d samtools-1.21.tar.bz2 && \
    tar -xvf samtools-1.21.tar && \
    cd samtools-1.21 && \
    ./configure --prefix=/software/samtools && \
    make -j 8 && \
    make install && \
    rm /software/samtools-1.21.tar && \
    conda create -n MAJIQ python=3.11 -y && \
    conda create -n SCSES python=3.11 -y
SHELL ["conda", "run", "-n", "SCSES", "/bin/bash", "-c"]
RUN pip install pandas numpy scipy scikit-learn && \
    pip install keras==2.15.0 && \
    pip install tensorflow==2.15.0.post1 && \
    pip install cython
    #conda install -c conda-forge gcc=12.1.0
# install MAJIQ in conda
SHELL ["conda", "run", "-n", "MAJIQ", "/bin/bash", "-c"]
RUN export HTSLIB_LIBRARY_DIR=/software/htslib/lib && \
    export HTSLIB_INCLUDE_DIR=/software/htslib/include && \
    chmod 777 -R /software/htslib && \
    cd /software && \
    conda install -c conda-forge gcc=12.1.0 && \
    #conda install conda-forge::zstd && \
    pip install git+https://bitbucket.org/biociphers/majiq_academic.git@v2.5.7 && \
    wget https://majiq.biociphers.org/app_download/majiq_license_academic_official.lic
# install STAR
SHELL ["/bin/bash", "-c"] 
RUN cd /software && \
    wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz -O Star.2.7.11.tar.gz && \
    tar -xzf Star.2.7.11.tar.gz && \
#   cd Star.2.7.11 && \
#   echo $PATH=/software/STAR-2.7.11b/bin/Linux_x86_64_static/STAR:$PATH >>/etc/profile && \
    rm Star.2.7.11.tar.gz
ENV PATH=$PATH:/software/STAR-2.7.11b/bin/Linux_x86_64_static/
SHELL ["conda", "run", "-n", "SCSES", "/bin/bash", "-c"]
RUN cd /software && \
    wget https://github.com/Xinglab/rmats-turbo/releases/download/v4.3.0/rmats_turbo_v4_3_0.tar.gz && \
    tar -zxvf rmats_turbo_v4_3_0.tar.gz && \
    cd rmats_turbo_v4_3_0 && \
    bash build_rmats && \
    rm /software/rmats_turbo_v4_3_0.tar.gz
    
ENV PATH=$PATH:$JAVA_HOME/bin:/software/rmats_turbo_v4_3_0/
# install MCR
RUN mkdir /MCR && \
    cd /MCR && \
    wget https://ssd.mathworks.com/supportfiles/downloads/R2022b/Release/10/deployment_files/installer/complete/glnxa64/MATLAB_Runtime_R2022b_Update_10_glnxa64.zip && \
    unzip -q MATLAB_Runtime_R2022b_Update_10_glnxa64.zip && \
    ./install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent && \
    cd / && \
    rm -rf MCR
# Install R packages
RUN R -e "install.packages(c('BiocManager','jsonlite','Matrix','reticulate','irlba','reshape2','R.matlab','hdf5r','R.oo','glmnet','caret','devtools'),dependencies=T,Ncpus=8)"
RUN R -e "devtools::install_version('Seurat', version = '4.4.0',upgrade='never')" 
#RUN R -e "devtools::install_github('dipterix/threeBrain',Ncpus=8,upgrade='never')"
#RUN R -e "devtools::install_github('jonclayden/RNifti',Ncpus=8,upgrade='never')" 
#RUN R -e "devtools::install_github('beauchamplab/raveio',Ncpus=8,upgrade='never')" 
RUN R -e "BiocManager::install(c('Rsamtools','Rhtslib','S4Vectors'),Ncpus=8,update=F)" 
RUN R -e "BiocManager::install(c('rtracklayer', 'Biostrings', 'GenomicRanges', 'IRanges', 'rhdf5'),Ncpus=8,update=F)" 
RUN R -e "devtools::install_version('BSgenome', version = '1.70.2', repos = 'https://bioconductor.org/packages/3.18/bioc',upgrade='never',dependencies=T)" 
RUN R -e "devtools::install_github('lvxuan12/SCSES',ref='SCSES_docker',Ncpus=8,upgrade='never')"
RUN R -e "install.packages('umap',dependencies=T,Ncpus=8)"
 #   echo PATH=/software/IRFinder-1.3.1/bin/util:$PATH >>/etc/profile
RUN cd /software && \
    wget https://github.com/RitchieLabIGH/IRFinder/archive/refs/tags/v2.0.1.tar.gz && \
    tar -zxvf v2.0.1.tar.gz
RUN mv ${CONDA_PREFIX}/lib/libstdc++.so.6 ${CONDA_PREFIX}/lib/libstdc++.so.6.bak &&\
    ln -s /usr/lib/x86_64-linux-gnu/libstdc++.so.6 ${CONDA_PREFIX}/lib/
ENV PATH=$PATH:/usr/lib/jvm/java-8-openjdk-amd64/bin:/software/IRFinder-2.0.1/bin/:/software/samtools/bin:/software/subread-2.0.6-source/bin
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/software/htslib/lib