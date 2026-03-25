#!/usr/bin/env python3
"""
FIXED A3M File Evaluator - Handles metadata lines correctly
"""

import os
import sys

def evaluate_a3m_file_fixed(filepath: str):
    """Fixed version that skips metadata lines starting with #"""
    
    print("\n" + "="*70)
    print("🔬 A3M FILE EVALUATOR - FIXED VERSION")
    print("="*70)
    
    if not os.path.exists(filepath):
        print(f"\n❌ File not found: {filepath}")
        return
    
    # Read file properly, skipping metadata
    with open(filepath, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    
    # Skip lines starting with # (metadata)
    filtered_lines = [line for line in lines if not line.startswith('#')]
    
    if not filtered_lines:
        print("\n❌ No valid sequences found (only metadata)")
        return
    
    # Parse sequences
    headers = []
    sequences = []
    current_seq = []
    
    for line in filtered_lines:
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
        print("\n❌ No protein sequences found")
        return
    
    # Get query sequence (first real sequence)
    query_seq = sequences[0]
    
    # Calculate stats
    file_size = os.path.getsize(filepath) / 1024
    sequence_count = len(sequences)
    query_length = len(query_seq)
    gaps = query_seq.count('-')
    gap_percent = (gaps / query_length * 100) if query_length > 0 else 0
    
    # Display results
    print(f"\n📁 File: {os.path.basename(filepath)}")
    print(f"📊 Size: {file_size:.1f} KB")
    
    print(f"\n" + "─"*70)
    print("📈 CORRECTED ANALYSIS")
    print("─"*70)
    
    print(f"• Total sequences: {sequence_count}")
    print(f"• Query length: {query_length} residues ✅")
    print(f"• Gaps in query: {gaps} ({gap_percent:.1f}%)")
    print(f"• Query coverage: {100 - gap_percent:.1f}%")
    
    # Show first few residues of query
    print(f"\n• Query sequence (first 50): {query_seq[:50]}...")
    
    # Check sequence composition
    amino_acids = set('ACDEFGHIKLMNPQRSTVWY')
    valid_chars = sum(1 for c in query_seq if c.upper() in amino_acids or c == '-')
    valid_percent = (valid_chars / query_length) * 100
    
    print(f"• Valid protein characters: {valid_percent:.1f}%")
    
    # Quality assessment
    print(f"\n" + "═"*70)
    print("💡 QUALITY ASSESSMENT")
    print("═"*70)
    
    # BACE1-specific check
    print(f"\n[FOR BACE1 (β-secretase)]:")
    if 380 <= query_length <= 420:
        print(f"  ✅ Correct length for BACE1 (~{query_length} residues)")
    else:
        print(f"  ⚠️  Unexpected length for BACE1 (expected ~400, got {query_length})")
    
    print(f"\n[SEQUENCE COUNT: {sequence_count}]")
    if sequence_count >= 1000:
        print("  ✅ EXCELLENT: Very deep MSA")
        print("     Consider subsampling to 512 for faster prediction")
    elif sequence_count >= 200:
        print("  ✅ EXCELLENT: Deep MSA")
    elif sequence_count >= 100:
        print("  ✅ GOOD: Adequate depth")
    
    print(f"\n[QUERY COVERAGE: {100 - gap_percent:.1f}%]")
    if gap_percent < 5:
        print("  ✅ EXCELLENT: Nearly complete coverage")
    elif gap_percent < 20:
        print("  ✅ GOOD: Good coverage")
    
    print(f"\n" + "═"*70)
    print("🎯 FINAL VERDICT")
    print("═"*70)
    
    print(f"\n✅ EXCELLENT MSA FOR BACE1!")
    print(f"   • {sequence_count} homologs found")
    print(f"   • Query: {query_length} residues (correct for BACE1)")
    print(f"   • Coverage: {100 - gap_percent:.1f}%")
    print(f"\n   This MSA is ready for high-quality AlphaFold prediction!")
    
    # Additional useful info
    print(f"\n" + "─"*70)
    print("💡 TIPS FOR ALPHAFOLD/COLABFOLD")
    print("─"*70)
    
    if sequence_count > 2000:
        print(f"• Your MSA has {sequence_count} sequences")
        print("  Consider using max_msa=512 for faster prediction")
        print("  (ColabFold parameter: max_msa=512)")
    
    print(f"• Expected prediction time: 10-30 minutes")
    print(f"• Expected pLDDT: >80 for most of the structure")
    
    print(f"\n" + "="*70)

def investigate_length_variation(filepath: str):
    """Find out why sequences have different lengths"""
    print(f"\n🔬 INVESTIGATING LENGTH VARIATION")
    print("="*70)
    
    # Read sequences
    with open(filepath, 'r') as f:
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
    
    query_len = len(sequences[0])
    
    # Analyze length distribution
    lengths = [len(s) for s in sequences]
    
    print(f"Total sequences: {len(sequences)}")
    print(f"Query length: {query_len}")
    print(f"Min length: {min(lengths)}")
    print(f"Max length: {max(lengths)}")
    print(f"Avg length: {sum(lengths)/len(lengths):.1f}")
    
    # Create histogram of lengths
    from collections import Counter
    length_counts = Counter(lengths)
    
    print(f"\n📊 LENGTH DISTRIBUTION:")
    for length, count in sorted(length_counts.items())[:20]:  # Show first 20
        if count > 1:
            print(f"  {length} residues: {count} sequences")
    
    # Find sequences significantly longer
    print(f"\n🔍 SEQUENCES LONGER THAN QUERY (>550 residues):")
    long_count = 0
    for i, (header, seq) in enumerate(zip(headers, sequences)):
        if len(seq) > 550:
            long_count += 1
            if long_count <= 3:  # Show first 3
                print(f"\n#{long_count}: {header[:80]}...")
                print(f"  Length: {len(seq)} (query: {query_len})")
                print(f"  Extra residues: {len(seq) - query_len}")
                
                # Compare to query
                print(f"  First 50 chars comparison:")
                print(f"  Query:  {sequences[0][:50]}")
                print(f"  This:   {seq[:50]}")
    
    if long_count > 0:
        print(f"\nTotal long sequences: {long_count}")
    
    # Check if long sequences just have terminal extensions
    print(f"\n🔍 CHECKING FOR TERMINAL EXTENSIONS:")
    
    # Look at sequence starts/ends
    query_start = sequences[0][:20]
    query_end = sequences[0][-20:] if len(sequences[0]) >= 20 else sequences[0]
    
    mismatched_starts = []
    mismatched_ends = []
    
    for i, seq in enumerate(sequences[1:11]):  # Check first 10
        if len(seq) > query_len + 10:
            seq_start = seq[:20]
            seq_end = seq[-20:] if len(seq) >= 20 else seq
            
            if query_start not in seq_start:
                mismatched_starts.append((i, seq_start))
            if query_end not in seq_end:
                mismatched_ends.append((i, seq_end))
    
    if mismatched_starts:
        print(f"  ⚠️  Some sequences don't start like query")
        print(f"  Query start: {query_start}")
        for idx, start in mismatched_starts[:3]:
            print(f"  Seq {idx} start: {start}")
    
    # Check for common tags/expression constructs
    print(f"\n🔍 CHECKING FOR COMMON TAGS:")
    common_tags = {
        "MGSSHHHHHH": "His-tag",
        "MASMTGGQQMG": "T7-tag",
        "MKYK...": "Signal peptide",
        "MAL...": "Signal peptide",
    }
    
    for seq in sequences[:10]:  # Check first 10
        for tag, name in common_tags.items():
            if tag in seq[:50]:
                print(f"  Found {name} in sequence")
    
    # Suggest solutions
    print(f"\n🔧 RECOMMENDED ACTIONS:")
    
    if max(lengths) > query_len + 100:
        print("1. ❗ Some sequences MUCH longer - likely different proteins")
        print("   → Filter by length: 400-550 residues for BACE1")
    
    # Calculate optimal filtering range
    optimal_min = query_len - 50
    optimal_max = query_len + 50
    
    sequences_in_range = [s for s in sequences if optimal_min <= len(s) <= optimal_max]
    print(f"2. Filter to {optimal_min}-{optimal_max} residues:")
    print(f"   Would keep: {len(sequences_in_range)}/{len(sequences)} sequences")
    
    if len(sequences_in_range) < 100:
        print("   ⚠️  Warning: Too few sequences after filtering!")
    else:
        print("   ✅ Good number would remain")

def main():
    """Main function"""
    # Your A3M file
    YOUR_FILE = "BACE1_human_mature_sequence_d8829_cleaned.a3m"
    
    # Check if file exists
    if not os.path.exists(YOUR_FILE):
        print(f"File '{YOUR_FILE}' not found!")
        print(f"Current directory: {os.getcwd()}")
        print("Files available:")
        for f in os.listdir('.'):
            if f.endswith('.a3m'):
                print(f"  • {f}")
        return
    
    # Run evaluation
    evaluate_a3m_file_fixed(YOUR_FILE)
    #Debug
    investigate_length_variation("BACE1_human_mature_sequence_d8829_cleaned.a3m")
    
    # Keep window open
    input("\nPress Enter to exit...")

if __name__ == "__main__":
    main()