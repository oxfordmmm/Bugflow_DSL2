#!/usr/bin/env bash

echo "building docker images matching *$1.Dockerfile"

base_image="docker/bugflow_base.Dockerfile"
echo "building base image $base_image"
docker build -f $base_image -t oxfordmmm/bugflow_base .
for file in $(find docker -iname "*$1.Dockerfile" -type f -maxdepth 1)
do
    if [ $file != $base_image ]
    then
        base_file=$(basename $file)
        image_name="${base_file%%.*}" 
        echo "Building oxfordmmm/${image_name} from $file"
        docker build -f $file -t oxfordmmm/${image_name} .
    else
        echo "not building $file, base image"
    fi
done
