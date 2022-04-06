#!/bin/bash

FQDIR="/mnt/ingolialab/kswain/NIKS015/NIKS016_extraseq/"
DATADIR="/mnt/ingolialab/ingolia/Tethering/NIKS018/"

CONSTANT="GATCCTGTAGCCCTAGACTTGATAGC"

BC_COUNT="/mnt/ingolialab/ingolia/Prog/barcode-assign/target/debug/bc-count"
BC_TABULATE="/mnt/ingolialab/ingolia/Prog/barcode-assign/bc-tabulate.py"

mkdir -p "${DATADIR}"

for FQ in `ls "${FQDIR}"/*.fastq | grep -v trim`
do
    SAMPLE=`basename "${FQ}" .fastq`

    OUTBASE="${DATADIR}/niks015-${SAMPLE}"
    COUNT="${OUTBASE}-count.txt"
    NBHD="${OUTBASE}-nbhd-count.txt"

    if [[ ! -e "${NBHD}" ]];
    then
        cutadapt -a "${CONSTANT}" \
	       --minimum-length 10 --discard-untrimmed \
	       "${FQ}" 2>"${DATADIR}/${SAMPLE}-cutadapt.txt" \
	  | "${BC_COUNT}" -f - -o "${COUNT}" -n "${OUTBASE}" &
        sleep 1
    else
        echo "${NBHD} exists, skipping ${SAMPLE}..."
    fi
done

wait

python "${BC_TABULATE}" --minsamples 2 -o "${DATADIR}/niks015.txt" `ls ${DATADIR}/niks015-*-nbhd-count.txt`
