#!/usr/bin/env python

"""
Standalone Python script to re-orient a circular GenBank file to start
at a specific gene (e.g., 'oriV', 'repA').
"""

import argparse
from Bio import SeqIO

def reorient_plasmid(record, gene_name):
    """
    Finds a gene and re-orients the circular plasmid record to start at that gene.
    Returns the new record, or the original record if the gene isn't found.
    """
    feature = None
    
    # Check for gene name in 'gene', 'product', or 'locus_tag' qualifiers
    for f in record.features:
        gene_names = f.qualifiers.get('gene', []) + \
                     f.qualifiers.get('product', []) + \
                     f.qualifiers.get('label', []) + \
                     f.qualifiers.get('locus_tag', [])
        
        if any(gene_name in name for name in gene_names):
            feature = f
            break
            
    if not feature:
        print(f"  WARNING: Gene '{gene_name}' not found in {record.id}. Saving original orientation.")
        return record # Return the original record if gene not found

    print(f"  Found '{gene_name}' at position {feature.location.start}. Re-orienting...")

    # If the gene is on the reverse strand, reverse-complement the whole plasmid
    # so the gene is on the positive strand at the start.
    if feature.location.strand == -1:
        print("  Gene is on reverse strand. Reverse-complementing plasmid.")
        
        # Store the original annotations before reverse-complementing
        original_annotations = record.annotations
        
        record = record.reverse_complement(id=True, name=True, description=True)
        
        # --- FIX: Manually copy the annotations to the new record ---
        record.annotations = original_annotations
        
        # We need to find the feature again in the new, reverse-complemented record
        for f in record.features:
            gene_names = f.qualifiers.get('gene', []) + \
                         f.qualifiers.get('product', []) + \
                         f.qualifiers.get('label', []) + \
                         f.qualifiers.get('locus_tag', [])
            if any(gene_name in name for name in gene_names):
                feature = f
                break
    
    # "Rotate" the circular plasmid so the gene is at the start
    start = feature.location.start
    
    # Ensure the new record maintains the original annotations and topology
    oriented_record = record[start:] + record[:start]
    oriented_record.id = record.id
    oriented_record.name = record.name
    oriented_record.description = record.description
    # This line now works because 'record' (whether original or rev-comped) has annotations
    oriented_record.annotations = record.annotations 
    
    return oriented_record

def main():
    """
    Main function to parse arguments and run the re-orientation.
    """
    # 1. Set up argument parser
    parser = argparse.ArgumentParser(description="Re-orient a circular GenBank plasmid to start at a specific gene.")
    parser.add_argument("-i", "--input", 
                        type=str, 
                        required=True, 
                        help="Input GenBank file (.gbk)")
    parser.add_argument("-o", "--output", 
                        type=str, 
                        required=True, 
                        help="Output file for the re-oriented GenBank record.")
    parser.add_argument("-g", "--gene", 
                        type=str, 
                        default="oriV",  # Set to 'oriV' as requested
                        help="Name of the gene to re-orient to (e.g., 'repA', 'oriV'). Default: 'oriV'")
    
    args = parser.parse_args()

    # 2. Read the input file
    print(f"Reading input file: {args.input}")
    try:
        record = SeqIO.read(args.input, "genbank")
    except FileNotFoundError:
        print(f"ERROR: Input file not found at {args.input}")
        return
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return

    # 3. Call the re-orient function
    print(f"Attempting to re-orient to gene: '{args.gene}'")
    oriented_record = reorient_plasmid(record, args.gene)

    # 4. Write the output file
    try:
        SeqIO.write(oriented_record, args.output, "genbank")
        print(f"\nSuccessfully wrote re-oriented file to: {args.output}")
    except Exception as e:
        print(f"An error occurred while writing the file: {e}")

if __name__ == "__main__":
    main()

