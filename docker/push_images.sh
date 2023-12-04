#!/usr/bin/env bash

echo "pushing docker images matching *$1.Dockerfile"

base_image="docker/bugflow_base.Dockerfile"
echo "Pushing base image $base_image"
docker push oxfordmmm/bugflow_base
for file in $(find docker -iname "*$1.Dockerfile" -type f -maxdepth 1)
do
    if [ $file != $base_image ]
    then
        base_file=$(basename $file)
        image_name="${base_file%%.*}" 
        echo "Pushing oxfordmmm/${image_name}"
        docker push oxfordmmm/${image_name}
    else
        echo "Not pushing $file, base image"
    fi
done
