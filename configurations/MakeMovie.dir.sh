FPS=15
EXT=png

IFS=$'\n'
DIR=`ls -d -Q --quoting-style=shell */`
DIR=`ls -d */`

for d in ${DIR};
do
  echo "->" $d
  mencoder "mf://$d/*.$EXT" -mf fps=$FPS -ovc lavc -lavcopts vcodec=msmpeg4v2:keyint=1:mbd=2:trell -o ${d/%?/.avi}
done

