#!/bin/bash

# Function to calculate proportion of variants
calculate_proportion() {
    vcf_file="$1"
    dp_threshold="${2:-10}"  # Set the default value of DP threshold to 10

    # Count total variants with FILTER equal to PASS
    total_pass_variants=$(zgrep -v '^#' "$vcf_file" | awk '$7 == "PASS"' | wc -l)

    # Count variants with FILTER equal to PASS and DP >= dp_threshold
    filtered_variants=$(zgrep -v '^#' "$vcf_file" | awk -F '\t' -v dp_thresh="$dp_threshold" '$7 == "PASS" {split($10, info, ":"); split(info[3], dp, ","); if(dp[1] >= dp_thresh) print $0}' | wc -l)

    # Calculate and display the proportion
    if [ "$total_pass_variants" -gt 0 ]; then
        proportion=$(echo "scale=4; $filtered_variants / $total_pass_variants" | bc)
        echo "Proportion of variants with FILTER=PASS and DP>=$dp_threshold: $proportion"
    else
        echo "No variants with FILTER=PASS found."
    fi
}

# Check if the VCF file is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <vcf_file.vcf.gz> [dp_threshold]"
    exit 1
fi

# Call the function with the provided VCF file and optional DP threshold
calculate_proportion "$1" "$2"
