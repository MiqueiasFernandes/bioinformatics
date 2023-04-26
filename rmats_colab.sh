#!/bin/bash
echo "*************instalando libs*************"
apt install wget samtools unzip curl python3 libgsl-dev -y 

echo "*************instalandh Hisat2****************"
wget -O hisat.zip https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download \
  && unzip -qq -o hisat.zip && mv hisat2-2.2.1/ hisat2

## instalar rmats
mkdir rmats && cd rmats && \
echo "*************instalando rMATS****************" && \
wget -O rmats "https://github.com/Xinglab/rmats-turbo/releases/download/v4.1.2/rmats_turbo_v4_1_2.tar.gz" && \
tar -xvf rmats && \
echo "*************compilando rMATS*************" && \
cd rmats_turbo* && make && \
cd .. && rm rmats && ln -s $(pwd)/$(ls rmats_turbo*/rmats.py) . 

echo "*************FINISHED*************"
