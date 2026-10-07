set -e
apt-get install -y -qq debootstrap
debootstrap --variant=minbase trixie /opt/trixie http://deb.debian.org/debian
chroot /opt/trixie bash -c 'apt-get update -qq && DEBIAN_FRONTEND=noninteractive apt-get install -y -qq ca-certificates gnupg wget'
chroot /opt/trixie bash -c 'mkdir -p /etc/apt/keyrings && wget -qO /etc/apt/keyrings/qgis-archive-keyring.gpg https://download.qgis.org/downloads/qgis-archive-keyring.gpg && printf "Types: deb\nURIs: https://qgis.org/debian\nSuites: trixie\nArchitectures: amd64\nComponents: main\nSigned-By: /etc/apt/keyrings/qgis-archive-keyring.gpg\n" > /etc/apt/sources.list.d/qgis.sources && apt-get update -qq && DEBIAN_FRONTEND=noninteractive apt-get install -y -qq qgis python3-qgis'
mkdir -p /opt/trixie/work
echo BUILD_DONE
