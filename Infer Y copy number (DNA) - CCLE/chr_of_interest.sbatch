#!/bin/bash

#SBATCH --job-name=chr_of_interest # Job name
#SBATCH -o slurm.main.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.main.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=ALL               # notifications 
#SBATCH --mail-user=nphaneuf@asu.edu # send-to address
#SBATCH -n 1
#SBATCH -p general
#SBATCH -q public
#SBATCH -t 0-23:00

cd /data/compres/CCLE_CURE_2023/DNA/chr_of_interest

module purge
module load fastqc-0.11.9-gcc-12.1.0

# load samtools
module load samtools-1.13-gcc-11.2.0

# Make directory if it does not exist
mkdir -p /scratch/nphaneuf/DNA_chr

# Specify the directory containing your BAM files
input_directory="/data/compres/CCLE_CURE_2023/DNA/chr_of_interest"

# Specify output directory
output_directory="/scratch/nphaneuf/DNA_chr"

# Assign CSV
csv_output="$output_directory/chr_of_interest.csv"

# Make directory if it does not exist
mkdir -p "$output_directory"

# Loop through files containing "chrY" in the filename
for bam_file in "$input_directory"/*chrX*.bam; do
    # Get the file name without the path
    filename=$(basename "$bam_file")

    # Run samtools cover and append to CSV and TXT, including the filename
    echo -n "$filename," >> "$csv_output"
    samtools coverage -r "X" "$filename" | tail -n +2 >> "$csv_output"  
 
done

for bam_file in "$input_directory"/*chrY*.bam; do
    # Get the file name without the path
    filename=$(basename "$bam_file")

    # Run samtools cover and append to CSV and TXT, including the filename
    echo -n "$filename," >> "$csv_output"
    samtools coverage -r "Y" "$filename" | tail -n +2 >> "$csv_output"  
 
done

for bam_file in "$input_directory"/*chr19*.bam; do
    # Get the file name without the path
    filename=$(basename "$bam_file")

    # Run samtools cover and append to CSV and TXT, including the filename
    echo -n "$filename," >> "$csv_output"
    samtools coverage -r "19" "$filename" | tail -n +2 >> "$csv_output"  
 
done

for bam_file in "$input_directory"/*chr8*.bam; do
    # Get the file name without the path
    filename=$(basename "$bam_file")

    # Run samtools cover and append to CSV and TXT, including the filename
    echo -n "$filename," >> "$csv_output"
    samtools coverage -r "8" "$filename" | tail -n +2 >> "$csv_output"  
 
done
