#!/usr/bin/env python3
"""
Script to replace 'start_time + <number>' patterns with 'start_time + n*dw' format.
Assumes dw = 200 (based on the codebase).
"""

import os
import re
from pathlib import Path

# Define the base dw value
DW = 200

def replace_time_offsets(content):
    """
    Replace patterns like 'start_time + 2*dw' with 'start_time + 2*dw'
    assuming dw = 200.
    """
    # Pattern to match 'start_time + <number>'
    pattern = r'start_time\s*\+\s*(\d+)'

    def replacement(match):
        offset = int(match.group(1))

        # Calculate the multiplier
        if offset % DW == 0:
            multiplier = offset // DW
            if multiplier == 1:
                return 'start_time + dw'
            else:
                return f'start_time + {multiplier}*dw'
        else:
            # If not divisible by dw, keep original
            print(f"Warning: {offset} is not divisible by {DW}, keeping original")
            return match.group(0)

    # Replace all occurrences
    new_content = re.sub(pattern, replacement, content)
    return new_content

def process_file(file_path):
    """Process a single file and replace patterns."""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()

        new_content = replace_time_offsets(content)

        # Only write if changes were made
        if new_content != content:
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(new_content)
            return True
        return False
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return False

def main():
    """Main function to process all Python files."""
    # Get the script's directory
    script_dir = Path(__file__).parent

    # File extensions to process
    extensions = ['.py', '.ipynb']

    # Track statistics
    files_modified = []
    files_scanned = 0

    print(f"Scanning directory: {script_dir}")
    print(f"Looking for patterns: start_time + <number>")
    print(f"Using dw = {DW}\n")

    # Walk through all files in the directory
    for root, dirs, files in os.walk(script_dir):
        # Skip hidden directories and common non-code directories
        dirs[:] = [d for d in dirs if not d.startswith('.') and d not in ['__pycache__', 'venv', 'env']]

        for file in files:
            file_path = Path(root) / file

            # Check if file has target extension
            if file_path.suffix in extensions:
                files_scanned += 1
                print(f"Processing: {file_path.relative_to(script_dir)}")

                if process_file(file_path):
                    files_modified.append(file_path.relative_to(script_dir))
                    print(f"  ✓ Modified")
                else:
                    print(f"  - No changes needed")

    # Print summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Files scanned: {files_scanned}")
    print(f"Files modified: {len(files_modified)}")

    if files_modified:
        print("\nModified files:")
        for file in files_modified:
            print(f"  - {file}")

    print("\nReplacement rules applied:")
    print(f"  start_time + {DW}   → start_time + dw")
    print(f"  start_time + {2*DW}  → start_time + 2*dw")
    print(f"  start_time + {3*DW}  → start_time + 3*dw")
    print(f"  start_time + {4*DW}  → start_time + 4*dw")
    print(f"  start_time + {5*DW}  → start_time + 5*dw")
    print("  etc...")

if __name__ == "__main__":
    main()
