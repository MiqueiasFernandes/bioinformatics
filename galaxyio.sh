#!/bin/bash

CRED=$1
FA=$2
FB=$3
KEY=`head -1 GalaxyEU_ids.txt`
USER=`tail -1 GalaxyEU_ids.txt`
API="https://usegalaxy.eu/api/datasets"

[ -z $FA ] && \
  echo "usage: galaxyio.sh GalaxyEU_ids.txt file.txt [galaxy_code]" \
  exit

[ -z $FB ] && echo "ENVIANDO $FA ..."
[ $FB ] && echo "RECEBENDO $FB ..."

 ## Enviar para o Galaxy.EU
[ -z $FB ] && \
  curl -T $FA --ssl ftp://ftp.usegalaxy.eu --user $USER

## Receber do Galaxy.EU
wget -O $FB "$API/$FA&key=$KEY"
