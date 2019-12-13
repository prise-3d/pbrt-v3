#! /bin/bash
if [ -z "$1" ]
  then
    echo "No argument supplied"
    echo "Need scenes folder path argument"
    exit 1
fi

if [ -z "$2" ]
  then
    echo "No argument supplied"
    echo "Need output folder path argument"
    exit 1
fi

numberofimages=10
numberofsamples=1

main_folder="$1/"
output_folder=$2
prefix="p3d_"

for folder in $(ls -d -- ${main_folder}*/)
do
  for file in $(ls $folder)
  do
    filename=$folder$file
    filename_fixed=${filename//\/\//\/}

    # check if filename contains 
    if [[ "$file" == ${prefix}* ]]; then
        echo ./pbrt --images ${numberofimages} --samples ${numberofsamples} --folder ${output_folder} ${filename_fixed}
    fi 
  done
done