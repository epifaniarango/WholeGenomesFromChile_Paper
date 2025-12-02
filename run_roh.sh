#!/usr/bin/env bash
#SBATCH --time=24:00:00
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1
#SBATCH --output=log/roh_%j.log

module load mamba

eval "$(conda shell.bash hook)"

# Activate SHAPEIT environment
conda activate bcftools


# Set input and output directories
input_dir="/scratch/earang/SMCJune/Population_VCFs"
output_dir="/scratch/earang/zhao2023/PopBoost/roh"



populations=("Yoruba" "French" "Han" "Huilliche-Chiloe" "Karitiana" "Lafkenche" "Pehuenche" "Quechua" "Surui")

for PopName in "${populations[@]}"; do
    bcftools roh -G30 --AF-dflt 0.4 "${input_dir}/${PopName}.vcf.gz" > "${output_dir}/${PopName}.txt"
done



combined_file="${output_dir}/combined_RG.txt"
echo -e "Sample\tChromosome\tStart\tEnd\tLength\tMarkers\tQuality\tPopulation" > $combined_file

for PopName in "${populations[@]}"; do
    grep "^RG" "${output_dir}/${PopName}.txt" | awk -v pop=$PopName '{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" pop}' >> $combined_file
done
