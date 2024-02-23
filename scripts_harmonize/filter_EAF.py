# filter the input file and exclude the rows that have no EAF value
import csv
import sys

def process_file(input, output_file_path):
    with open(input, 'r', newline='', encoding='utf-8') as file:
        count = 0
        header = file.readline().strip().split('\t')

        eaf_index = header.index('EAF') if 'EAF' in header else -1 # find EAF column index

        if eaf_index != -1:
            with open(output_file_path, 'w', newline='', encoding='utf-8') as output:
                tsv_writer = csv.writer(output, delimiter='\t')
                tsv_writer.writerow(header)

                for line in file:
                        columns = line.strip().split('\t')
                        # check if column exists and EAF is not NA or .
                        if len(columns) > eaf_index:
                        
                            if columns[eaf_index] != "NA" and columns[eaf_index] != '.':
                                output.write(line) # only the valid snps are included in the output file
                            else:
                                count += 1


    print("Number of excluded Snps because of unvalid EAF: ", count)
    print("Results are safed in", output_file_path)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("parameters are missing")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    process_file(input_file, output_file) 