

cp /etc/apt/sources.list /etc/apt/sources.list.bkp
echo -e "
####################################################################################\n\
#                               Repositórios Oficiais                              #\n\
####################################################################################\n\
\n\
deb http://deb.debian.org/debian/ stable main contrib non-free\n\
# deb-src http://deb.debian.org/debian/ stable main contrib non-free\n\
\n\
deb http://deb.debian.org/debian/ stable-updates main contrib non-free\n\
# deb-src http://deb.debian.org/debian/ stable-updates main contrib non-free\n\
\n\
deb http://deb.debian.org/debian-security stable/updates main contrib non-free\n\
# deb-src http://deb.debian.org/debian-security stable/updates main contrib non-free\n\
\n\
## Debian Stretch Backports\n\
# deb http://ftp.debian.org/debian stretch-backports main contrib non-free\n\
# deb-src http://ftp.debian.org/debian stretch-backports main contrib non-free\n\
\n\
####################################################################################\n\n\
deb http://cran.rstudio.com/bin/linux/debian stretch-cran34/\n\
\n\
" > /etc/apt/sources.list


cp /etc/network/interfaces /etc/network/interfaces.bkp
echo -e "\n\
auto eno1\n\
iface eno1 inet static\n\
address 192.168.38.207\n\
netmask 255.255.255.0\n\
broadcast 192.168.38.255\n\
network 192.168.38.1\n\
gatway 192.168.38.1\n\
\n\
auto eno2\n\
iface eno2 inet static\n\
address 169.0.0.20\n\
netmask 255.255.255.0\n\
\n\
auto enp3s0f1\n\
iface enp3s0f1 inet static\n\
address 169.0.0.21\n\
netmask 255.255.255.0\n\
\n\
auto enp3s0f0\n\
iface enp3s0f0 inet static\n\
address 169.0.0.22\n\
netmask 255.255.255.0\n" > /etc/network/interfaces


cp /etc/hosts /etc/hosts.bkp
echo -e "\n\n\
169.0.0.10 bioserver1\n\
169.0.0.20 bioserver2\n\
169.0.0.30 bioserver3\n\
169.0.0.40 bioserver4\n\
\n" >> /etc/hosts




apt update
apt upgrade
apt install firmware-linux firmware-linux-nonfree locales synaptic gdebi arc arj cabextract lhasa p7zip p7zip-full p7zip-rar rar unrar unace unzip xz-utils zip default-jre default-jdk aptitude build-essential most cups-pdf poppler-utils python3-pip build-essential libssl-dev libffi-dev python-dev python3-venv dirmngr aptitude texlive-xetex resolvconf ncbi-blast+ libText-Soundex.perl git repeatmasker-recon libJSON.perl libURI.perl liblwp-useragent-determined-perl screen augustus bamtools samtools zlib1g-dev libbamtools-dev libboost-iostreams-dev libboost-all-dev libhts-dev bcftools tabix libsam-dev libbam-dev libghc-bzlib-dev liblzma-dev libssl-dev libcurl4-openssl-dev genometools htop iftop apache2 curl openssh-server openssh-client mysql-server php libapache2-mod-php php-mcrypt php-mysql 


apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF'
apt install r-base

## en_US ISO-8859-1
## en_US.UTF-8 UTF-8
## pt_BR ISO-8859-1
## pt_BR.UTF-8 UTF-8
## DEFAULT: pt_BR.UTF-8
dpkg-reconfigure locales

echo "domain ufes.br\nsearch ufes.br\n"> /etc/resolvconf/resolv.conf.d/head

echo "nameserver 192.168.56.16" > /etc/resolvconf/resolv.conf.d/tail
echo "nameserver 192.168.56.17" >> /etc/resolvconf/resolv.conf.d/tail
dpkg-reconfigure resolvconf


pip3 install pandas biopython jupyter biocode gffutils bcbio-gff

####https://peteris.rocks/blog/sun-grid-engine-installation-on-ubuntu-server/

##apt install gridengine-master gridengine-client

##apt install gridengine-exec




update-rc.d apache2 disable
systemctl set-default multi-user.target
echo -e "\n\nalias gui=\"su - -c 'service gdm3 start'\"\n\n" >> /home/cluster/.bashrc


###configure sshfs
#rm ~/.ssh/ -r
ssh-keygen -b 4096
ssh-copy-id root@bioserver1
mkdir /home/cluster/shared && chmod 777 /home/cluster/shared
sshfs root@bioserver1:/home/cluster/shared /home/cluster/shared -o allow_other




