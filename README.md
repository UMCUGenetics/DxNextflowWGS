# DxNextflowWGS [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4709755.svg)](https://doi.org/10.5281/zenodo.4709755)

Genome Diagnostics Nextflow WGS workflow

## Get Nextflow Modules and install OpenJDK and Nextflow
```bash
sh install.sh
```

## Running WGS workflow
```bash
nextflow run WGS.nf -c WGS.config --fastq_path <fastq_dir_path> --outdir <output_dir_path> --email <email> [-profile slurm|mac]
```

## QDNAseq GUIX container
```bash
guixr pack -f squashfs -RR -S /bin=bin r r-qdnaseq-hmf r-getoptlong perl bash glibc-utf8-locales tzdata coreutils procps grep sed bootstrap-binaries
cp /guix/path.gz.squashfs QDNAseq_v1.9.2-HMF.1.gz.squashfs
unsquashfs QDNAseq_v1.9.2-HMF.1.gz.squashfs
chmod -R 0775 squashfs-root
singularity build QDNAseq_v1.9.2-HMF.1.sif squashfs-root
```

## BAF container
```bash
guixr pack -f squashfs -RR -S /bin=bin bio-vcf r r-ggplot2 r-gtools r-pastecs bash glibc-utf8-locales tzdata coreutils procps grep sed bootstrap-binaries
cp /guix/path.gz.squashfs baf_nextflow.gz.squashfs
unsquashfs baf_nextflow.gz.squashfs
chmod -R 0775 squashfs-root
singularity build baf_nextflow.sif squashfs-root
```
