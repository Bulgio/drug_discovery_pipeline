"""
CLEAN and evaluate A3M file - essential for Boltzmann generators
"""

import os
import sys

def clean_a3m_file(input_path: str, output_path: str = None):
    """
    Remove metadata lines and fix A3M format
    Returns path to cleaned file
    """
    if output_path is None:
        output_path = input_path.replace('.a3m', '_cleaned.a3m')
    
    print(f"\n🔧 Cleaning A3M file: {os.path.basename(input_path)}")
    
    with open(input_path, 'r') as f:
        lines = f.readlines()

    # Filter out problematic lines
    cleaned_lines = []
    for line in lines:
        line = line.rstrip('\n')

        # Skip metadata lines
        if line.startswith('#'):
            print(f"   Removing metadata: {line[:50]}...")
            continue
        
        # Skip empty lines
        if not line.strip():
            continue
        
        # Keep everything else
        cleaned_lines.append(line)
    
    # Write cleaned file
    with open(output_path, 'w') as f:
        for line in cleaned_lines:
            f.write(line + '\n')
    
    print(f"   Original: {len(lines)} lines")
    print(f"   Cleaned: {len(cleaned_lines)} lines")
    print(f"   Saved to: {output_path}")
    
    return output_path

def filter_long_seq(input_path: str,n = 50):
    """
    Filtering sequences too long to make biological sense, 50+- residues than the query
    """


    outclean_path = input_path.replace('.a3m', f'_{n}_filtered.a3m')
    outliers_path = input_path.replace('.a3m', f'_{n}_outliers.a3m')
    
    print(f"\n🔧 Filtering A3M file: {os.path.basename(input_path)}")

        # Read sequences
    with open(input_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    

    sequences = []
    headers = []
    current_seq = []
    
    for line in lines:
        if line.startswith('>'):
            if current_seq:
                sequences.append(''.join(current_seq))
                current_seq = []
            headers.append(line)
        else:
            current_seq.append(line)
    
    if current_seq:
        sequences.append(''.join(current_seq))
    
    if len(sequences) < 2:
        return
    #Lenght of the query sequence
    query_len = len(sequences[0])

    outliers = []
    cleaned = []

    for header,sequence in zip(headers,sequences):
        seq_len=len(sequence)
        if not (query_len - n) <= seq_len <= (query_len + n):
            # Put sequences in a separate container
            outliers.append(header)
            outliers.append(sequence)
        else:
            cleaned.append(header)
            cleaned.append(sequence)

    with open(outclean_path, 'w') as f:
        for item in cleaned:
            f.write(item + '\n')
        print(f'Filtered file saved in {outclean_path}')

    with open(outliers_path, 'w') as f:
        for item in outliers:
            f.write(item + '\n')
        print(f'Outliers saved in {outliers_path}')

        # Print summary
    total_sequences = len(sequences)
    kept_sequences = len(cleaned) // 2  # Each sequence has header + sequence
    removed_sequences = len(outliers) // 2
    
    print(f"\n📊 FILTERING RESULTS:")
    print(f"  Total sequences: {total_sequences}")
    print(f"  Kept: {kept_sequences} sequences")
    print(f"  Removed: {removed_sequences} sequences")

    return outclean_path, outliers_path

def validate_a3m_for_boltzmann(filepath: str):
    """
    Validate A3M file is suitable for Boltzmann generators
    """
    print(f"\n" + "="*70)
    print("🧪 VALIDATION FOR BOLTZMANN GENERATORS")
    print("="*70)
    
    # Read and parse
    with open(filepath, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    
    headers = []
    sequences = []
    current_seq = []
    
    for line in lines:
        if line.startswith('>'):
            if current_seq:
                sequences.append(''.join(current_seq))
                current_seq = []
            headers.append(line)
        else:
            current_seq.append(line)
    
    if current_seq:
        sequences.append(''.join(current_seq))
    
    if not sequences:
        print("❌ ERROR: No sequences found after cleaning!")
        return False
    
    query_seq = sequences[0]
    query_len = len(query_seq)
    
    print(f"\n📊 VALIDATION CHECKS:")
    
    # Check 1: Query sequence sanity
    print(f"\n1. Query sequence check:")
    print(f"   Length: {query_len} residues")
    
    # Count valid amino acids
    aa_set = set('ACDEFGHIKLMNPQRSTVWY')
    valid_chars = sum(1 for c in query_seq if c.upper() in aa_set or c == '-')
    invalid_chars = query_len - valid_chars
    
    if invalid_chars == 0:
        print(f"   ✅ All characters are valid AAs or gaps")
    else:
        print(f"   ⚠️  {invalid_chars} invalid characters found")
        invalid = [c for c in query_seq if c.upper() not in aa_set and c != '-']
        print(f"   Invalid chars: {set(invalid)}")
    
    # Check 2: Sequence length consistency
    print(f"\n2. Length consistency check:")
    lengths = [len(seq) for seq in sequences]
    if all(l == query_len for l in lengths):
        print(f"   ✅ All {len(sequences)} sequences have same length")
    else:
        print(f"   ⚠️  Lengths vary: {min(lengths)} to {max(lengths)}")
        print(f"   This may cause issues with some Boltzmann implementations")
    
    # Check 3: Gap patterns
    print(f"\n3. Gap analysis:")
    query_gaps = query_seq.count('-')
    gap_percent = (query_gaps / query_len) * 100
    
    # Check if gaps are reasonable
    if gap_percent < 10:
        print(f"   ✅ Low gap percentage: {gap_percent:.1f}%")
    elif gap_percent < 30:
        print(f"   ⚠️  Moderate gaps: {gap_percent:.1f}%")
    else:
        print(f"   ❌ High gaps: {gap_percent:.1f}% - may affect modeling")
    
    # Check 4: For Boltzmann specific needs
    print(f"\n4. Boltzmann generator considerations:")
    
    # Deep MSAs are good for Boltzmann
    if len(sequences) >= 1000:
        print(f"   ✅ Deep MSA ({len(sequences)} seqs) - good for diversity")
    elif len(sequences) >= 100:
        print(f"   ✅ Adequate depth ({len(sequences)} seqs)")
    else:
        print(f"   ⚠️  Shallow MSA ({len(sequences)} seqs) - limited diversity")
    
    # Check sequence diversity
    unique_seqs = len(set(sequences))
    diversity = unique_seqs / len(sequences)
    print(f"   Sequence diversity: {diversity:.2f}")
    
    # Check 5: Format compliance
    print(f"\n5. Format compliance:")
    problematic_headers = [h for h in headers if ' ' in h and len(h.split()) > 2]
    if problematic_headers:
        print(f"   ⚠️  Complex headers found (may need simplification)")
        print(f"   Example: {problematic_headers[0][:50]}...")
    else:
        print(f"   ✅ Headers look clean")
    
    print(f"\n" + "="*70)
    print("🎯 RECOMMENDATION FOR BOLTZMANN:")
    print("="*70)
    
    if query_len > 50 and len(sequences) >= 100 and gap_percent < 30:
        print(f"\n✅ SUITABLE for Boltzmann generators")
        print(f"   • Query: {query_len} residues")
        print(f"   • Sequences: {len(sequences)} homologs")
        print(f"   • Gaps: {gap_percent:.1f}%")
    else:
        print(f"\n⚠️  MAY NEED PREPROCESSING")
        if gap_percent > 30:
            print(f"   • Consider removing columns with >50% gaps")
        if len(sequences) < 100:
            print(f"   • Consider deeper MSA search")
    
    return True

def main():
    """Clean and validate A3M for Boltzmann use"""
    
    # Your file
    input_file = "BACE1_human_mature_sequence_d8829.a3m"
    
    if not os.path.exists(input_file):
        print(f"File not found: {input_file}")
        return
    
    # Step 1: Clean the file
    cleaned_file = clean_a3m_file(input_file)

    # Step 2: Filter the file
    filtered_file,outliers_file=filter_long_seq(cleaned_file,25)
    
    # Step 2: Validate for Boltzmann
    validate_a3m_for_boltzmann(filtered_file)
    
    # Step 3: Show what was fixed
    print(f"\n" + "="*70)
    print("📝 WHAT WAS FIXED:")
    print("="*70)
    print(f"""
    
    Files created:
    • {input_file} - Original (keep as backup)
    • {cleaned_file} - Cleaned version
    • {filtered_file} - Use this 
    • {outliers_file} - Outliers sequences and headers, keep this for documentation
    """)
    
    # Keep window open
    input("\nPress Enter to exit...")

if __name__ == "__main__":
    main()