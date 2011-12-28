#!/bin/bash

#cluster-fork mkdir /tmp/temp
DIR=/tmp/temp
FILE=$1
thishost=`hostname --short`

if [ "$thishost" == "levlabhn" ]; then
  scp $DIR/$FILE compute-0-0:$DIR
  ssh -f compute-0-0 FastCopy.sh $FILE
  scp $DIR/$FILE compute-0-5:$DIR
  scp $DIR/$FILE compute-0-6:$DIR
  scp $DIR/$FILE compute-0-7:$DIR
  scp $DIR/$FILE compute-0-8:$DIR
fi

if [ "$thishost" == "compute-0-0" ]; then
  scp $DIR/$FILE compute-0-1:$DIR
  scp $DIR/$FILE compute-0-2:$DIR
  scp $DIR/$FILE compute-0-3:$DIR
  scp $DIR/$FILE compute-0-4:$DIR
  wall "DONE: $0 $@"
fi
