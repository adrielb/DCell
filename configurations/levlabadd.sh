#!/bin/bash

export USERNAME=$1

useradd ${USERNAME}

export DBNAME=db_${USERNAME}
export DBUSER=db_${USERNAME}
export DBPASS=`./pwdgen.sh 12`

mysql << EOF
CREATE DATABASE ${DBNAME};
GRANT ALL PRIVILEGES ON ${DBNAME}.* TO ${DBUSER}@localhost IDENTIFIED BY "${DBPASS}";
FLUSH PRIVILEGES;
EXIT
EOF

su ${USERNAME}

echo -e "Database Name: ${DBNAME}\nDatabase user: ${DBUSER}\nPassword: ${DBPASS}\n" >> ~/mysql-access
chmod 600 ~/mysql-access

add cron job for mysqlhotcopy -u ${DBUSER} -p ${DBPASS} ${DBNAME} ~/db-backup

mkdir ~/.ssh
chmod 700 ~/.ssh
ssh-keygen -q -f ~/.ssh/id_rsa -t rsa -N ''

mkdir ~/public_html
mkdir ~/public_html/public_data

wget -P /tmp/ http://wordpress.org/latest.tar.gz
tar zxvf latest.tar.gz
mv wordpress/* public_html/
rmdir wordpress

sed -e "\
s/putyourdbnamehere/${DBNAME}/;\
s/usernamehere/${DBUSER}/;\
s/yourpasswordhere/${DBPASS}/" \
~/public_html/wp-config-sample.php > ~/public_html/wp-config.php

wget -P /tmp/ https://api.wordpress.org/secret-key/1.1/
KEYS=`cat /tmp/index.html`
rm /tmp/index.html

chgrp apache /home/abergman/public_html/wp-content
chmod 775 /home/abergman/public_html/wp-content


echo "Now visit http://128.220.136.155/~${USERNAME}/wp-admin/install.php"

pip3pi3k
