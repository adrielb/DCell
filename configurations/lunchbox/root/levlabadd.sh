#!/bin/bash

web=/home/$1/public_html
dat=${web}/public_data

useradd $1
chmod 711 /home/$1
mkdir     ${web}
chmod a+w ${web}
chown $1:$1 ${web}
mkdir     ${dat}
chown $1:$1 ${dat}
ln -s ${dat} /home/$1/
passwd $1
