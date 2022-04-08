#!/bin/bash

DATADIR="./work"

CONSTANT="GATCCTGTAGCCCTAGACTTGATAGC"

BARCODE_ASSIGN="./barcode-assign/target/debug"
BC_COUNT="${BARCODE_ASSIGN}/bc-count"
BC_TABULATE="${BARCODE_ASSIGN}/bc-tabulate"

mkdir -p "${DATADIR}"

for SRR in `cut -f1 samples.txt`
do
    SAMPLE=`grep "${SRR}" samples.txt | cut -f2`

    OUTBASE="${DATADIR}/${SAMPLE}"
    COUNT="${OUTBASE}-count.txt"
    NBHD="${OUTBASE}-count.txt"

    if [[ ! -e "${NBHD}" ]];
    then
        fastq-dump -Z "${SRR}" 2>"${DATADIR}/${SAMPLE}-fastqdump.txt" \
	  | cutadapt -a "${CONSTANT}" \
		   --minimum-length 10 --discard-untrimmed \
		   - 2>"${DATADIR}/${SAMPLE}-cutadapt.txt" \
	  | "${BC_COUNT}" -f - -o "${COUNT}" -n "${OUTBASE}" &
        sleep 1
    else
        echo "${NBHD} exists, skipping ${SAMPLE}..."
    fi
done

wait

"${BC_TABULATE}" --minsamples 2 -o "${DATADIR}/niks015.txt" `ls ${DATADIR}/niks015-*-count.txt`
"${BC_TABULATE}" --minsamples 2 -o "${DATADIR}/niks018.txt" `ls ${DATADIR}/niks018-*-count.txt`
