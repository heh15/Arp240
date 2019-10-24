#!/bin/bash

work_dir=/1/home/heh15/workingspace/Arp240/
image_dir=$work_dir'ratio/NGC5257/startImage'
filenames=$work_dir'*CO*/finalImage/NGC5257*combine_contsub.image'

for filename in $filenames
do
cp -r $filename $image_dir
done


