#!/bin/bash
set -euo pipefail

workflow_path='/hpc/diaggen/software/production/DxNextflowWGS'

# Set input and output dirs
input=`realpath -e $1`
output=`realpath $2`
email=$3
mkdir -p $output && cd $output
mkdir -p log

if ! { [ -f 'workflow.running' ] || [ -f 'workflow.done' ] || [ -f 'workflow.failed' ]; }; then
touch workflow.running

sbatch <<EOT
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --mem 5G
#SBATCH --gres=tmpspace:10G
#SBATCH --job-name Nextflow_WGS
#SBATCH -o log/slurm_nextflow_wgs.%j.out
#SBATCH -e log/slurm_nextflow_wgs.%j.err
#SBATCH --mail-user $email
#SBATCH --mail-type FAIL
#SBATCH --export=NONE
#SBATCH --account=diaggen

module load Java/1.8.0_60

/hpc/diaggen/software/tools/nextflow run $workflow_path/WGS.nf \
-c $workflow_path/WES.config \
--fastq_path $input \
--outdir $output \
--email $email \
-profile slurm \
-resume -ansi-log false

if [ \$? -eq 0 ]; then
    echo "Nextflow done."

    echo "WGS workflow completed successfully."
    rm workflow.running
    touch workflow.done

    exit 0
else
    echo "Nextflow failed"
    rm workflow.running
    touch workflow.failed
    exit 1
fi
EOT
else
echo "Workflow job not submitted, please check $output for 'workflow.status' files."
fi