#!/bin/bash

DATADIR=" ~/Tethering/PacBio190731"
CCS="${DATADIR}/m64044_190724_081616.ccs.bam"
GENOME="${DATADIR}/saccharomyces_cerevisiae.fa"

BARCODE_ASSIGN="./barcode-assign/target/debug"

"${BARCODE_ASSIGN}/bc-pbr" \
    -i "${CCS}" \
    -l lib-specs.txt \
    -o "${DATADIR}/pacbio-190731"

/mnt/ingolialab/linux-x86_64/src/smrtlink_7.0.1.66975/smrtcmds/bin/blasr \
    --hitPolicy allbest --bam --nproc 36 \
    "${DATADIR}/pacbio-190731-frags.fasta" \
    "${GENOME}" \
    --out "${DATADIR}/pacbio-190731-frags-aligned.bam"

"${BARCODE_ASSIGN}/bc-pbj" \
    -i "${DATADIR}/pacbio-190731" \
    -o "${DATADIR}/pacbio-190731"
