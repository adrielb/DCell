#!/bin/bash
# request Bourne again shell as shell for job
#$ -S /bin/bash
# merge stderr with stdout
#$ -j y
# assign a name to submitted job
#$ -N MovieMaking
# change output file name
#$ -o output/$JOB_NAME.$JOB_ID

if [ "${HOSTNAME%%-*}" != "compute" ]
then
  echo "Do not run on head node"
  exit 1
fi

REMOTE_PATH=$1
REMOTE_DIR=${REMOTE_PATH%/*}
REMOTE_FILE=${REMOTE_PATH##/*/}
MOVIE_FILE=${REMOTE_FILE/%spl/avi}

BITRATE=5000
FPS=30

echo "Remote path: $REMOTE_PATH"
echo "Remote dir:  $REMOTE_DIR"
echo "Remote file: $REMOTE_FILE"
echo "Movie file:  $MOVIE_FILE"

#Copy spool file to local
#SpoolToBytes to break spool into gray
#Convert gray into jpg
#Assemble mpeg from jpg
#Copy mpeg file to precis

. ~/.bashrc

set -x
BASH_IFS=$IFS

cd $TMP

REMOTE_DIR="\"$REMOTE_DIR\""
IFS=""

scp -i ~/PrecisID/id_rsa_precis LevLabHN@10.16.107.22:\"$REMOTE_PATH\" $TMP && \
SpoolToBytes "$REMOTE_FILE" && \
cd DIR."$REMOTE_FILE" && \
for f in *.gray; do RawToJPEG $f ${f/%gray/jpg}; done && \
mencoder "mf://*.jpg" -quiet -mf fps=$FPS -o ../"$MOVIE_FILE" -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=$BITRATE && \
cd .. && \
scp -i ~/PrecisID/id_rsa_precis "$MOVIE_FILE" LevLabHN@10.16.107.22:$REMOTE_DIR


#using ImageMagick 'convert'
#for f in *.gray; do convert -size 512x512 -depth 16 -define quantum:format=unsigned-integer -quality 100 $f $f.jpg; done

#removing file extension
#find . -name "*.gray" | cut -f1-4 -d'.'
