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


nthreads=$(less /proc/cpuinfo | grep processor | wc -l)

if [ -z "$3" ]
  then
    echo "No nthreads argument supplied, by default use ${nthreads} cpu ressources.."
else
  nthreads=$3
fi

numberofimages=1
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
        ./pbrt --images ${numberofimages} --samples ${numberofsamples} --folder ${output_folder} --nthreads ${nthreads} --normals 1 --zbuffer 1 ${filename_fixed}
    fi 
  done
done
