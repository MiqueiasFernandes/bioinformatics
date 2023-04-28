#!/bin/bash

### Sample file => GalaxyEU_ids.txt
#   email@host.net
#   2348342bcf23465cb233323cc64235d234a
### end sample file => GalaxyEU_ids.txt

CRED=$1
FA=$2
FB=$3
KEY=`head -1 $CRED`
USER=`tail -1 $CRED`
API="https://usegalaxy.eu/api/datasets"

[ -z $FA ] && \
  echo "usage: galaxyio.sh GalaxyEU_ids.txt file.txt [galaxy_code]" && \
  exit

[ -z $FB ] && echo "ENVIANDO $FA ..."
[ $FB ] && echo "RECEBENDO $FB ..."


## Enviar para o Galaxy.EU
[ -z $FB ] && \
  curl -T $FA --ssl ftp://ftp.usegalaxy.eu --user $USER

## Receber do Galaxy.EU
[ $FB ] && wget -O $FB "$API/$FA/display?key=$KEY"

