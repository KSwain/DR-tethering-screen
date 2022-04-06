#! /bin/bash

FQGZDIR="/mnt/ingolialab/FastQ/190515_50SR_HS4KA/Swain/"
DATADIR="/mnt/ingolialab/ingolia/Tethering/NIKS018/"

CONSTANT="GATCCTGTAGCCCTAGACTTGATAGC"

BC_COUNT="/mnt/ingolialab/ingolia/Prog/barcode-assign/target/debug/bc-count"
BC_TABULATE="/mnt/ingolialab/ingolia/Prog/barcode-assign/bc-tabulate.py"

mkdir -p "${DATADIR}"

for FQGZ in `ls "${FQGZDIR}"/*.fastq.gz`
do
    SAMPLE=`basename "${FQGZ}" | sed s/_S._L005_R1_001.fastq.gz//`

    OUTBASE="${DATADIR}/niks018-${SAMPLE}"
    COUNT="${OUTBASE}-count.txt"
    NBHD="${OUTBASE}-nbhd-count.txt"

    if [[ ! -e "${NBHD}" ]];
    then
        cutadapt -a "${CONSTANT}" \
	       --minimum-length 10 --discard-untrimmed \
	       "${FQGZ}" 2>"${DATADIR}/${SAMPLE}-cutadapt.txt" \
	  | "${BC_COUNT}" -f - -o "${COUNT}" -n "${OUTBASE}" &
        sleep 1
    else
        echo "${NBHD} exists, skipping ${SAMPLE}..."
    fi
done

wait

python "${BC_TABULATE}" --minsamples 2 -o "${DATADIR}/niks018.txt" `ls ${DATADIR}/niks018-*-nbhd-count.txt`
