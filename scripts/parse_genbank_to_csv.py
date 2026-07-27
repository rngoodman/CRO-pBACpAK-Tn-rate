#!/usr/bin/env python

"""
Standalone Python script to parse one or more GenBank files and
extract plasmid information and feature annotations into two
separate CSV files for use in R/circlize.
"""

import argparse
import os
import csv
from Bio import SeqIO

def main():
    """
    Main function to parse arguments and run the GenBank parsing.
    """
    # 1. Set up argument parser
    parser = argparse.ArgumentParser(
        description="Parse GenBank files to create CSVs for circlize."
    )
    parser.add_argument(
        "-i", "--input",
        type=str,
        nargs='+',  # Accept one or more input files
        required=True,
        help="One or more input GenBank files (.gbk)."
    )
    parser.add_argument(
        "-o", "--output_prefix",
        type=str,
        required=True,
        help="Output prefix for the CSV files (e.g., 'my_data')."
    )
    
    args = parser.parse_args()

    # Lists to hold all data
    plasmid_data = []
    gene_data = []
    
    print(f"\nExtracting data from {len(args.input)} file(s)...")

    # 2. Loop through each input file and parse data
    for filepath in args.input:
        if not os.path.exists(filepath):
            print(f"  WARNING: File not found: {filepath}. Skipping.")
            continue
        
        try:
            # --- NEW: Get plasmid name from filename ---
            base_name = os.path.basename(filepath)
            plasmid_name = os.path.splitext(base_name)[0]
            
            record = SeqIO.read(filepath, "genbank")
            print(f"  Processing: {plasmid_name} (from file {filepath})")

            # 1. Add plasmid size to plasmid_df
            plasmid_data.append({
                'name': plasmid_name, # Use filename-derived name
                'start': 0,
                'end': len(record.seq)
            })
            
            # 2. Add gene features to gene_df
            for feature in record.features:
                # You can add more feature types here if needed
                if feature.type in ('CDS', 'gene', 'rep_origin', 'misc_feature'):
                    
                    # --- Get gene name ---
                    gene_name = feature.qualifiers.get('gene', 
                                  feature.qualifiers.get('label', 
                                      feature.qualifiers.get('locus_tag', 
                                          feature.qualifiers.get('product', ['unknown'])
                                      )
                                  )
                              )[0]
                    
                    # --- Get feature type (with override) ---
                    feature_type = feature.type  # Get the default type
                    label_value = feature.qualifiers.get('label', [None])[0]
                    
                    if label_value in ("λ repressor", "-lambda repressor (fragment)", "-lambdarepressor", "lambdarepressor", "λ repressor (fragment)", "-lambda repressor", "lambda repressor"):
                        feature_type = "repressor"

                    if gene_name in ("tnpA", "tnpR", "AmpR promoter", "AmpR (fragment)", "insL1", "ISKpn25", "IRR", "IRL", "IS10R", "Tn1000", "TNPA_ECOLI", "TNR1_ECOLI", "ALTA2_ALTAL", "TNPX_ECOLI (fragment)", "IS5", "IS1A", "Tn2"):
                        feature_type = "insert"

                    if gene_name in ("mcr1", "mcr-1"):
                        feature_type = "mcr1"

                    if gene_name in ("KanR", "KanR"):
                        feature_type = "aph(3′)-Ia"

                    gene_data.append({
                        'plasmid': plasmid_name, # Use filename-derived name
                        'start': int(feature.location.start),
                        'end': int(feature.location.end),
                        'direction': '=>' if feature.location.strand == 1 else '<=' if feature.location.strand == -1 else '==',
                        'gene_name': gene_name,
                        'type': feature_type, # Use the (potentially overridden) type
                        'product': feature.qualifiers.get('product', ['unknown'])[0]
                    })
        except Exception as e:
            print(f"  ERROR processing file {filepath}: {e}")

    # 3. Define output filenames
    plasmids_csv_path = f"{args.output_prefix}_plasmids.csv"
    genes_csv_path = f"{args.output_prefix}_genes.csv"

    # 4. Write CSV files
    try:
        # Write plasmids_df.csv
        with open(plasmids_csv_path, 'w', newline='') as f:
            # --- CHANGED: Renamed 'name' to 'plasmid' ---
            writer = csv.DictWriter(f, fieldnames=['name', 'start', 'end'])
            writer.writeheader()
            writer.writerows(plasmid_data)
        
        # Write genes_df.csv
        with open(genes_csv_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['plasmid', 'start', 'end', 'direction', 'gene_name', 'type', 'product'])
            writer.writeheader()
            writer.writerows(gene_data)
            
        print("\nDone.")
        print(f"Created: {plasmids_csv_path}")
        print(f"Created: {genes_csv_path}")

    except Exception as e:
        print(f"\nAn error occurred while writing the CSV files: {e}")

if __name__ == "__main__":
    main()

