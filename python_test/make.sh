#!/bin/bash
python skrypt.py
time parallel ../build/src/Raytracing -x 1920 -y 1080 ::: output/*
# time parallel -j 1515 ../build/src/Raytracing -x 1920 -y 1080 ::: output/*
# ffmpeg -r 50 -f image2 -s 1920x1080 -i ./output/%04d.bmp -vcodec libx264 -crf 25 -pix_fmt yuv420p test.mp4
# cp -r output/*.bmp /mnt/c/Linux/BMP/
rm -rf output/*
# mv test.mp4 /mnt/c/Linux/