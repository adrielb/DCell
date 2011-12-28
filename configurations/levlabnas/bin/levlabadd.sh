#!/usr/local/bin/bash

if [ "$1" == "" ]; then exit 0; fi

zfs create tank/home/$1
chown $1:levlab /mnt/tank/home/$1
chmod 700 /mnt/tank/home/$1


