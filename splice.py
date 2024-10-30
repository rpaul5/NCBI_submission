import csv
import os

def read_fasta(file):
    """Generator function to read a FASTA file."""
    header = ""
    sequence = []
    
    with open(file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    yield header, "".join(sequence)
                header = line
                sequence = []
            else:
                sequence.append(line)
        yield header, "".join(sequence)

def splice_contigs(csv_file):
    # Load the splice information from the CSV, grouping by old_file, sequence_name, and new_file
    contigs_to_splice = {}
    with open(csv_file, 'r') as csvfile:
        reader = csv.reader(csvfile)
        header_row = next(reader)  # Skip the header row
        for row in reader:
            try:
                old_file = row[0].strip()
                sequence_name = row[1].strip()
                start = int(row[2].strip())
                end = int(row[3].strip())
                new_file = row[4].strip()
                
                if not os.path.exists(old_file):
                    print(f"Error: {old_file} does not exist.")
                    continue
                
                if old_file not in contigs_to_splice:
                    contigs_to_splice[old_file] = {}
                if new_file not in contigs_to_splice[old_file]:
                    contigs_to_splice[old_file][new_file] = {}
                if sequence_name not in contigs_to_splice[old_file][new_file]:
                    contigs_to_splice[old_file][new_file][sequence_name] = []
                
                contigs_to_splice[old_file][new_file][sequence_name].append((start, end))
            except ValueError as e:
                print(f"Error processing row: {row}. Error: {e}")
                continue
    
    # Process each old_file and apply all splices, writing the result to each specified new_file
    for old_file, new_files in contigs_to_splice.items():
        for new_file, sequences in new_files.items():
            with open(new_file, 'w') as out:
                print(f"Processing {old_file}, outputting to {new_file}...")
                for header, sequence in read_fasta(old_file):
                    sequence_name = header[1:]  # Get sequence name without '>'
                    
                    if sequence_name in sequences:
                        # Sort regions in reverse order to avoid index shifting
                        splice_regions = sequences[sequence_name]
                        splice_regions.sort(reverse=True, key=lambda x: x[0])
                        
                        spliced_sequence = sequence
                        for start, end in splice_regions:
                            spliced_sequence = spliced_sequence[:start-1] + spliced_sequence[end:]
                            print(f"Splicing {sequence_name} in {old_file} from {start} to {end} for output in {new_file}")
                        
                        out.write(f"{header}\n")
                        out.write(f"{spliced_sequence}\n")
                    else:
                        # Write unmodified sequence if no splicing is specified for this sequence_name
                        out.write(f"{header}\n")
                        out.write(f"{sequence}\n")
                        print(f"No splicing for {sequence_name} in {old_file}; copied to {new_file}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python splice_contigs.py regions.csv")
        sys.exit(1)
    
    csv_file = sys.argv[1]
    splice_contigs(csv_file)

