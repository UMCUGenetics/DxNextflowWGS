#!/bin/bash
set -euo pipefail

workflow_path='/hpc/diaggen/software/production/DxNextflowWGS'

# Set input and output dirs
input=`realpath -e $1`
output=`realpath $2`
email=$3
optional_params=( "${@:4}" )
mkdir -p $output && cd $output
mkdir -p log

if ! { [ -f 'workflow.running' ] || [ -f 'workflow.done' ] || [ -f 'workflow.failed' ]; }; then
touch workflow.running

file="${output}/log/nextflow_trace.txt"

if [ -e "${file}" ]; then
    # Extract the current suffix from the file name
    current_suffix=$(echo "${file}" | grep -oE '[0-9]+$')

    if [ -z "${current_suffix}" ]; then
        # If no suffix found, set it to 0
        current_suffix=0
    fi

    # Increment the suffix
    new_suffix=$((current_suffix + 1))

    # Create the new file name with the incremented suffix
    new_file="${file%.*}_$new_suffix.${file##*.}"

    # Rename the file
    mv "${file}" "${new_file}"
fi

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

export NXF_JAVA_HOME='/hpc/diaggen/software/tools/jdk-20.0.2/'

/hpc/diaggen/software/tools/nextflow run $workflow_path/WGS.nf \
-c $workflow_path/WGS.config \
--fastq_path $input \
--outdir $output \
--email $email \
-profile slurm \
-resume -ansi-log false \
${optional_params[@]:-""}

if [ \$? -eq 0 ]; then
    echo "Nextflow done."

    echo "Zip work directory"
    find work -type f | egrep "\.(command|exitcode)" | zip -@ -q work.zip

    echo "Remove work directory"
    rm -r work

    echo "Creating md5sum"
    find -type f -not -iname 'md5sum.txt' -exec md5sum {} \; > md5sum.txt

    echo "WGS workflow completed successfully."
    rm workflow.running
    touch workflow.done

    echo "Change permissions"
    chmod 775 -R $output
    
    exit 0
else
    echo "Nextflow failed"
    rm workflow.running
    touch workflow.failed
    
    echo "Change permissions"
    chmod 775 -R $output
    
    exit 1
fi
EOT
else
echo "Workflow job not submitted, please check $output for 'workflow.status' files."
fi
