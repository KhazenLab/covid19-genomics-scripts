#!/bin/bash

date=$(date +%d%m%Y)
ref=EPI_ISL_402125
INDS=($(for i in SubFiles/${date}/*_${date}.fasta; do echo $(basename ${i%.*}); done))

mkdir AnalysisFiles
mkdir AnalysisFiles/${date}

for qry in ${INDS[@]}; do

    echo ${qry}
    nucmer -p nucmer_${ref}_${qry} Reference/${ref}.fasta SubFiles/${date}/${qry}.fasta
    show-coords -r -c -l nucmer_${ref}_${qry}.delta > AnalysisFiles/${date}/nucmer_${ref}_${qry}.coords
    show-snps -C nucmer_${ref}_${qry}.delta > AnalysisFiles/${date}/nucmer_${ref}_${qry}.snps
    dnadiff  -d nucmer_${ref}_${qry}.delta -p AnalysisFiles/${date}/dnadiff_nucmer_${ref}_${qry}

    mv *.delta AnalysisFiles/${date}
    mv *.mgaps AnalysisFiles/${date}
    mv *.ntref AnalysisFiles/${date}


done

zip -r AnalysisFiles/${date}.zip AnalysisFiles/${date}



