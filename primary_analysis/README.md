# Barcode-to-fragment assignment

Required files:
* circular consensus sequencing (CCS) file `m64044_190724_081616.ccs.bam`
* yeast genome Fasta file `saccharomyces_cerevisiae.fa`
* [`barcode-assign` tools](https://github.com/ingolia-lab/barcode-assign)

## `pacbio-barcodes.sh`

Generates barcode-to-fragment files

From `bc-pbr`:
* `pacbio-190731-frags.fasta`
* `pacbio-190731-read-inserts-good.txt`
* `pacbio-190731-read-fates.txt`

From `bc-pbj`:
* `pacbio-190731-read-aligns-all.txt`
* `pacbio-190731-read-aligns-unique.txt`
* `pacbio-190731-barcode-assign-all.txt`
* `pacbio-190731-barcode-assign-umabig.txt`
* `pacbio-190731-barcode-assign-unique.txt`
* `pacbio-190731-barcode-assign.bed`

## `pacbio-assign.sh` and `pacbio-assign.R`

Shell script generates:
* `pacbio-190731-facs-assign.bed`
* `pacbio-190731-facs-assign-gene.bed`
and then runs the R script to calculate some library statistics and generate:
* `pacbio-190731-facs-assign-yorf.txt`

# Barcode counting

Required files:
* [`NCBI SRA Toolkit`](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software) to download FastQ data
* [`barcode-assign` tools](https://github.com/ingolia-lab/barcode-assign)
* [`cutadapt`](https://cutadapt.readthedocs.io/en/stable/index.html)

## `count-barcodes.sh`

Processes raw sequencing data to count barcodes:
* `niks015.txt` and `niks018.txt`

## `frag-counts.R`

Uses `pacbio-190731-barcode-assign-all.txt` and
`pacbio-190731-facs-assign-yorf.txt` to associate a fragment with the
barcodes in the count tables:
* `niks018-assigned-counts.csv` and `niks015-assigned-counts.csv`


