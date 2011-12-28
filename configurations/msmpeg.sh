#!/bin/bash

#input
PREFIX=img
SUFFIX=jpg
FPS=5
OUTPUT=output

#set -e

rm -rf $OUTPUT.avi divx2pass.log

img=$PREFIX???.$SUFFIX

width=`identify $img | awk '{print $3}' | awk -Fx '{print $1}'`
height=`identify $img | awk '{print $3}' | awk -Fx '{print $2}'`

echo $width $height

#
# compute the optimal bitrate 
#optimal_bitrate = block_bit_rate * 25 * width * height / 256
#
bbr=60
#obr=`expr $width \* $height \* $bbr \* 25 / 256`

#OPTS="vbitrate=2000:mbd=2:trell:dia=2:keyint=1 -vf scale=-1 -msglevel all=-1"
OPTS="vbitrate=4000"

#IMAGELIST="mf://$PREFIX.?.$SUFFIX mf://$PREFIX.??.$SUFFIX mf://$PREFIX.???.$SUFFIX"
IMAGELIST="mf://$PREFIX???.$SUFFIX"

echo "Pass 1"
mencoder $IMAGELIST -ovc lavc -lavcopts vcodec=msmpeg4v2:vpass=1:$OPTS -mf type=$SUFFIX:w=$width:h=$height:fps=$FPS -nosound -o /dev/null
echo "Pass 2"
mencoder $IMAGELIST -ovc lavc -lavcopts vcodec=msmpeg4v2:vpass=2:$OPTS -mf type=$SUFFIX:w=$width:h=$height:fps=$FPS -nosound -o $OUTPUT.avi

#mencoder $IMAGELIST -ovc lavc -lavcopts vcodec=msmpeg4v2:$OPTS -mf type=$SUFFIX:w=$width:h=$height:fps=$FPS -nosound -o $OUTPUT.avi

ls -hl $OUTPUT.avi
