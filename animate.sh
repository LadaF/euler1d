#!/bin/bash


for i in *.dat ; do
 sed "s/calc1/`echo $i | sed "s/\.dat$//"`/g" anim.plot | gnuplot
done
 
ffmpeg -r 10 -i %03d.jpg anim.mp4
exit 0
