#!/bin/bash

# goal of this script is to look inside a folder
# and extract pertinent data from all of its
# subdirectories


echo "Output of file in f1.txt"
exec 1> f1.txt

path1="/home/kazuumi/lus/project/chem761_final"
path2="/home/ruisun/proj/ruisun/B0"
PATH="$PATH:$path1:$path2"

grep 'w lines' $path2/3/fort.8
