#!/bin/bash

while [ 0 ]
do 
	rsync -a /home/abergman/dicty/figs/* lunchbox:public_data/figs/
	date
	sleep 1
done

