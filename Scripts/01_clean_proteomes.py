import os
import re
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def clean_protein_fasta(input_file, output_file, summary_data):
    """
    Cleans a protein FASTA file for use with OrthoFinder and tracks statistics.
    """
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY*")  # Standard 20 amino acids + stop codon
    cleaned_records = []
    
    # Initialize counters for this file
    initial_count = 0
    removed_invalid_chars = 0
    removed_short_seq = 0
    removed_high_gap = 0
    
    for record in SeqIO.parse(input_file, "fasta"):
        initial_count += 1
        seq_str = str(record.seq).upper()
        
        # Track original sequence for comparison
        original_length = len(seq_str)
        original_gap_count = seq_str.count('-') + seq_str.count('X')
        
        # 1. Clean the sequence: uppercase, replace ambiguous 'X' with gap
        seq_str = seq_str.replace('X', '-')
        
        # Count and remove invalid amino acid characters
        valid_seq_chars = []
        invalid_char_count = 0
        for aa in seq_str:
            if aa in valid_aa:
                valid_seq_chars.append(aa)
            else:
                invalid_char_count += 1
                valid_seq_chars.append('-')  # Replace invalid with gap
        
        if invalid_char_count > 0:
            removed_invalid_chars += 1
        
        seq_str = ''.join(valid_seq_chars)
        
        # 2. Filter sequences
        # Calculate gap percentage
        gap_count = seq_str.count('-')
        gap_percentage = gap_count / len(seq_str) if len(seq_str) > 0 else 1.0
        
        # Skip sequences that are too short
        if len(seq_str) < 20:
            removed_short_seq += 1
            continue
        
        # Skip sequences with too many gaps
        if gap_percentage > 0.5:
            removed_high_gap += 1
            continue
        
        # 3. Clean the header
        header = record.description
        new_id = header.split()[0]
        new_id = re.sub(r'[|:,;()\[\]]', '_', new_id)
        
        # 4. Create cleaned record
        cleaned_record = SeqRecord(
            Seq(seq_str),
            id=new_id,
            description=""
        )
        cleaned_records.append(cleaned_record)
    
    # Write cleaned sequences to file
    final_count = len(cleaned_records)
    if cleaned_records:
        SeqIO.write(cleaned_records, output_file, "fasta")
    
    # Add to summary data
    summary_data.append({
        'Filename': os.path.basename(input_file),
        'Initial_Genes': initial_count,
        'Removed_Invalid_Chars': removed_invalid_chars,
        'Removed_Short_Seq': removed_short_seq,
        'Removed_High_Gap': removed_high_gap,
        'Final_Genes': final_count,
        'Genes_Lost': initial_count - final_count,
        'Retention_Rate (%)': round((final_count / initial_count * 100), 2) if initial_count > 0 else 0
    })
    
    return final_count

def batch_clean_fasta(directory):
    """
    Cleans all FASTA files in a given directory and generates a summary CSV.
    """
    summary_data = []
    
    for filename in os.listdir(directory):
        if filename.endswith(('.faa', '.fasta', '.fa', '.pep')):
            input_path = os.path.join(directory, filename)
            base_name = os.path.splitext(filename)[0]
            output_path = os.path.join(directory, f"cleaned_{base_name}.faa")
            
            print(f"Processing {filename}...")
            final_count = clean_protein_fasta(input_path, output_path, summary_data)
            print(f"  Initial: {summary_data[-1]['Initial_Genes']}, Final: {final_count}")
    
    # Save summary to CSV
    csv_path = os.path.join(directory, "cleaning_summary.csv")
    with open(csv_path, 'w', newline='') as csvfile:
        if summary_data:
            fieldnames = summary_data[0].keys()
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(summary_data)
    
    print(f"\nSummary saved to: {csv_path}")
    
    # Print summary table
    print("\n" + "="*80)
    print("CLEANING SUMMARY")
    print("="*80)
    print(f"{'Filename':<25} {'Initial':<8} {'Final':<8} {'Lost':<8} {'Retention %':<12}")
    print("-"*80)
    
    total_initial = 0
    total_final = 0
    
    for item in summary_data:
        print(f"{item['Filename'][:24]:<25} {item['Initial_Genes']:<8} {item['Final_Genes']:<8} "
              f"{item['Genes_Lost']:<8} {item['Retention_Rate (%)']:<12}")
        total_initial += item['Initial_Genes']
        total_final += item['Final_Genes']
    
    print("-"*80)
    total_lost = total_initial - total_final
    total_retention = round((total_final / total_initial * 100), 2) if total_initial > 0 else 0
    print(f"{'TOTAL':<25} {total_initial:<8} {total_final:<8} {total_lost:<8} {total_retention:<12}")
    print("="*80)

# --- Run the script ---
if __name__ == "__main__":
    target_directory = r"E:\fish\orth"
    
    # Check if directory exists
    if not os.path.exists(target_directory):
        print(f"Error: Directory '{target_directory}' does not exist.")
    else:
        print(f"Starting batch cleaning in: {target_directory}")
        batch_clean_fasta(target_directory)
        print("\nBatch cleaning complete!")