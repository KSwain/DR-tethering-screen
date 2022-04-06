#!/bin/bash

DATADIR=~/Tethering/PacBio190731

ASSIGN_BED="${DATADIR}/pacbio-190731-barcode-assign.bed"
FACS_BED="${DATADIR}/pacbio-190731-facs-assign.bed"

GENOME_BED=sac_cer_yassour_clean.bed

grep '_FACS' "${ASSIGN_BED}" \
    | sed 's/_FACS//' \
	> "${DATADIR}/pacbio-190731-facs-assign.bed"
bedtools intersect -s -wa -wb \
         -a "${FACS_BED}" -b "${GENOME_BED}" \
         > "${DATADIR}/pacbio-190731-facs-assign-gene.bed"

R --no-save < pacbio-assign.R
