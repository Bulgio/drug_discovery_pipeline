#!/usr/bin/env python3
"""
Fix PDB format for RDKit compatibility.
Handles:
1. 4-character residue names (e.g., LIG1 → LIG)
2. Element symbols in wrong columns (moves from col 79-80 to col 77-78)
Ensures each line is exactly 80 characters with element symbols in columns 77-78.
"""
import sys

def fix_pdb_for_rdkit(input_file, output_file):
    print(f"Fixing PDB file: {input_file} -> {output_file}")
    print("=" * 60)

    with open(input_file, 'r') as f:
        lines = f.readlines()

    fixed_lines = []
    atom_counter = 1
    fixed_residue_count = 0
    fixed_element_count = 0

    for line_num, line in enumerate(lines, 1):
        original_line = line.rstrip('\n\r')
        
        if original_line.startswith(('HETATM', 'ATOM')):
            print(f"\n--- Processing atom {atom_counter} (line {line_num}) ---")
            print(f"Original: '{original_line}'")
            print(f"Length: {len(original_line)}")
            
            # Step 1: Fix residue name (LIG1 → LIG)
            # Check if there's a 4-char residue name at cols 17-21
            if len(original_line) > 21:
                resname_raw = original_line[17:21]
                if resname_raw and len(resname_raw.strip()) > 3:
                    # Found 4-char residue name
                    old_res = resname_raw
                    new_res = resname_raw[:3]  # Take first 3 chars
                    print(f"  Found 4-char residue: '{old_res}' -> '{new_res}'")
                    
                    # Rebuild line with 3-char residue
                    # Format: take up to col 16, add new residue (3 chars), then rest from col 21
                    line = original_line[:17] + new_res + original_line[21:]
                    fixed_residue_count += 1
                else:
                    line = original_line
            else:
                line = original_line
            
            # Step 2: Fix element position (move from col 79-80 to col 77-78)
            if len(line) >= 80:
                # Element is at cols 79-80 (index 78-79)
                element = line[78:80].strip()
                print(f"  Element found at cols 79-80: '{element}'")
                
                # Take first 76 chars
                prefix = line[:76]
                
                # Format element (ensure 2 chars, right-justified)
                if len(element) == 1:
                    element = f" {element}"
                elif len(element) > 2:
                    element = element[:2]
                elif len(element) == 0:
                    # Try to guess from atom name
                    atom_name = line[12:16].strip()
                    element = ''.join(c for c in atom_name if c.isalpha())[:2]
                    if len(element) == 1:
                        element = f" {element}"
                    print(f"  Guessed element from atom name: '{element.strip()}'")
                
                # New line: prefix (76) + element (2)
                new_line = prefix + element
                fixed_element_count += 1
                print(f"  Element moved to cols 77-78: '{new_line[76:78]}'")
            else:
                # Line too short, pad it
                print(f"  Line too short ({len(line)}), padding to 76 chars")
                new_line = line
                while len(new_line) < 76:
                    new_line += " "
                
                # Add default element
                atom_name = line[12:16].strip() if len(line) > 12 else ""
                element = ''.join(c for c in atom_name if c.isalpha())[:2]
                if len(element) == 1:
                    element = f" {element}"
                elif len(element) == 0:
                    element = " C"
                
                new_line = new_line + element
                print(f"  Added element '{element.strip()}' at cols 77-78")
            
            # Ensure exactly 80 characters
            if len(new_line) > 80:
                print(f"  Truncating from {len(new_line)} to 80")
                new_line = new_line[:80]
            elif len(new_line) < 80:
                print(f"  Padding from {len(new_line)} to 80")
                new_line = new_line + " " * (80 - len(new_line))
            
            print(f"Final line ({len(new_line)} chars):")
            print(f"  '{new_line}'")
            print(f"  Element at cols 77-78: '{new_line[76:78]}'")
            
            fixed_lines.append(new_line)
            atom_counter += 1

        elif line.startswith('END'):
            fixed_lines.append('END')
        else:
            fixed_lines.append(original_line)

    # Write fixed file
    print("\n" + "=" * 60)
    print(f"Writing fixed file: {output_file}")
    
    with open(output_file, 'w') as f:
        for i, line in enumerate(fixed_lines):
            if i < len(fixed_lines) - 1:
                f.write(line + '\n')
            else:
                f.write(line)

    print(f"\nSummary:")
    print(f"  Atoms processed: {atom_counter-1}")
    print(f"  Residue names fixed: {fixed_residue_count}")
    print(f"  Element positions fixed: {fixed_element_count}")
    
    # Verification
    print("\n" + "=" * 60)
    print("Verification (first 3 atoms):")
    with open(output_file, 'r') as f:
        for i, line in enumerate(f.readlines()[:3]):
            line = line.rstrip('\n')
            print(f"\nAtom {i+1}:")
            print(f"  Length: {len(line)}")
            if len(line) >= 78:
                print(f"  Residue (cols 18-20): '{line[17:20]}'")
                print(f"  Element (cols 77-78): '{line[76:78]}'")
                print(f"  Full line: '{line}'")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: fix_rdkit.py <input.pdb> <output.pdb>")
        print("Example: fix_rdkit.py ligand_raw.pdb ligand_fixed.pdb")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    fix_pdb_for_rdkit(input_file, output_file)