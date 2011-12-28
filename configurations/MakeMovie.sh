FPS=30
EXT=png

mencoder "mf://*.$EXT" -mf fps=$FPS -ovc lavc -lavcopts vcodec=msmpeg4v2:keyint=1:mbd=2:trell -o movie.avi

