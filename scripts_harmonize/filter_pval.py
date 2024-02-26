import gzip
import sys

# alternative to filtering with awk

def process_file(file_path, output_file):
    count = 0
    try:
       with open(file_path, 'r') as file:
            header = file.readline().strip().split('\t') 

            # Find Index for Pvalue
            column_pval_index = header.index('P-value') if 'P-value' in header else -1 # if needed adapt to column name

            if column_pval_index != -1: # index is -1 if column not found
                with open(output_file, 'w') as output:
                    for line in file:
                        columns = line.strip().split('\t')

                        # check if column exists and pvalue is <= 0.05 & not NA
                        if len(columns) > column_pval_index and columns[column_pval_index] != "NA":
                            if float(columns[column_pval_index]) <= 0.05:
                                count += 1
                                output.write(line)

    except FileNotFoundError:
        print(f'File {file_path} not found.')

    except Exception as e:
        print(f'Exception occured: {e}')

    print(f'Associations with pval <= 0.05: {count}')

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("parameters are missing")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    process_file(input_file, output_file)