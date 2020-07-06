#!/bin/sh
date=19062020

java -Xmx4g -jar snpEff/snpEff.jar NC_045512.2 mutations_${date}.vcf > mutations_annotations_${date}.vcf

