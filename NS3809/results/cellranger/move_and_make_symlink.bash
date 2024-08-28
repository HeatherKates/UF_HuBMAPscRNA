#!/bin/bash

for dir in NS-*; do
    mv "$dir" /orange/timgarrett/Campbell_Thompson/scRNA/cellranger_output/
    ln -s /orange/timgarrett/Campbell_Thompson/scRNA/cellranger_output/"$dir" "$dir"
done
