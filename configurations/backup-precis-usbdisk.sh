#!/bin/bash

SOURCE=$HOME
DEST=/media/usbdisk/
EXCLUDE=$HOME/Research/configurations/rsync.exclude
LOG=$HOME/backup.log

START_TIME=`date +%s`

rsync --verbose --stats --archive --exclude-from=$EXCLUDE $SOURCE $DEST > $LOG 2>&1

END_TIME=`date +%s`

echo "TIME:" `expr $END_TIME - $START_TIME` "sec" >> $LOG
