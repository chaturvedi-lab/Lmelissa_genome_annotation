import argparse

# Define the range (10,000 bp) for matching positions
range_threshold = 10000

# Create a dictionary to store functional annotations
annotations = {}

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Match SNP data with functional annotations and save the results to an output file.")
parser.add_argument("snp_file", help="Path to the SNP file")
parser.add_argument("annotation_file", help="Path to the functional annotation file")
parser.add_argument("output_file", help="Path to the output file")
args = parser.parse_args()

# Read the functional annotation file
with open(args.annotation_file, 'r') as annotation_file:
    for line in annotation_file:
        parts = line.strip().split('\t')
        # Split the first column before "HRSCAF" and take the scaffold part
        scaffold = parts[0].split('_HRSCAF_')[0]
        start, stop = int(parts[1]), int(parts[2])
        # Take the last four columns
        annotation_data = parts[-4:]
        go_info = parts[6] if len(parts) > 6 else ""  # Information from the 7th column if available
        transposase_info = parts[8] if len(parts) > 8 else ""  # Information from the 9th column if available
        annotations[scaffold] = (start, stop, annotation_data, go_info, transposase_info)

# Process the SNP file and save the results to the output file
with open(args.output_file, 'w') as output_file:
    with open(args.snp_file, 'r') as snp_file:
        for line in snp_file:
            parts = line.strip().split('\t')
            scaffold, position = parts[0], int(parts[1])

            if scaffold in annotations:
                start, stop, annotation_data, go_info, transposase_info = annotations[scaffold]
                if (start <= position <= stop) or (position - start <= range_threshold) or (stop - position <= range_threshold):
                    go_numbers = list(set([col for col in annotation_data if "GO" in col]))
                    ipr_numbers = [col for col in annotation_data if "IPR" in col]

                    matched_data = '\t'.join(parts + [",".join(go_numbers), ",".join(ipr_numbers), go_info, transposase_info])
                    output_file.write(matched_data + '\n')
                else:
                    # No match found; you can handle this case as needed
                    output_file.write('\t'.join(parts + ['No Match']) + '\n')
            else:
                # No match found; you can handle this case as needed
                output_file.write('\t'.join(parts + ['No Match']) + '\n')

