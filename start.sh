#!/bin/sh

echo "configurando eno1 em dhcp"

dhclient eno1


echo "inicializando servico sshfs"


sshfs root@bioserver1:/home/cluster/shared /home/cluster/shared -o allow_other -f


echo "serviço sshfs terminado com sucesso"

exit 0
