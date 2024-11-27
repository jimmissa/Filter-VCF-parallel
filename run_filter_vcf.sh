#!/bin/bash

# Input variables
VCF_PATH='../../../species_wide_vcf_v6/data/species_wide_vcf_v6/species_wide.v6.vcf.gz'
OUT_FILE="filtered_output.vcf"
SNP=1
BIALLELIC=0
INDEL=2
DP=3
GQ=20
AD=0
COVG_FILE="average_depth_new.txt"
FRAC_CALLED=0.5
QUANT=0.9
CHROMOSOMES="chromosomes.txt"
SAMPLES=("BMS-1" "CB6-3" "CC1009" "CC1010" "CC1373" "CC1952" "CC2342" "CC2343" "CC2344" "CC2931" "CC2932" "CC2935" "CC2936" "CC2937" "CC2938" "CC3059" "CC3060" "CC3061" "CC3062" "CC3063" "CC3065" "CC3068" "CC3069" "CC3071" "CC3073" "CC3076" "CC3079" "CC3084" "CC3086" "CC3268" "CL3-1" "GB117" "GB119" "GB123" "GB13" "GB138" "GB141" "GB66" "HH93-01" "JG4-1" "MW46-01" "MW46-02" "MW46-03" "MW46-04" "MW46-05" "MW46-06" "MW46-08" "MW46-09" "MW46-10" "MW46-12" "MW46-16" "MW46-17" "MW46-19" "MW46-20" "MW47-02" "PW10-01" "PW10-03" "PW10-06" "PW21-01" "UTM7-1" "W13-1" "W13-2")
PROCESSES=-1

# Check required files
if [[ ! -f "$VCF_PATH" ]]; then
    echo "Error: VCF file not found at $VCF_PATH."
    exit 1
fi

if [[ ! -f "$COVG_FILE" ]]; then
    echo "Error: Coverage file not found at $COVG_FILE."
    exit 1
fi

if [[ ! -f "$CHROMOSOMES" ]]; then
    echo "Error: Chromosomes file not found at $CHROMOSOMES."
    exit 1
fi

# Execute Python script
python3 filter_vcf.py \
    --out "$OUT_FILE" \
    --vcf_path "$VCF_PATH" \
    --snp "$SNP" \
    --biallelic "$BIALLELIC" \
    --indel "$INDEL" \
    --DP "$DP" \
    --GQ "$GQ" \
    --AD "$AD" \
    --covg_file "$COVG_FILE" \
    --frac_called "$FRAC_CALLED" \
    --quant "$QUANT" \
    --samples "${SAMPLES[@]}" \
    --chromosomes "$CHROMOSOMES" \
    --processes "$PROCESSES"
