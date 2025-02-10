import os
import sys
import gzip
import pandas as pd

def read_barcodes(file_path):
    df = pd.read_csv(file_path)
    return dict(zip(df['barcode'], df['cell_type']))

def split_fastq(fastq_dir, sample, output_dir, barcodes_csv, chunk_size=100000):
    barcodes = read_barcodes(barcodes_csv)
    forward_fq = f"{fastq_dir}/{sample}_1.final.trimmed.fastq.gz" # найти способ выбирать нужную пару файлов в заданной папке по префиксу
    reverse_fq = f"{fastq_dir}/{sample}_2.final.trimmed.fastq.gz"
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    forward_handles = {cell_type: gzip.open(f"{output_dir}/{sample}_{cell_type}_1.fastq.gz", "wt") for cell_type in set(barcodes.values())}
    reverse_handles = {cell_type: gzip.open(f"{output_dir}/{sample}_{cell_type}_2.fastq.gz", "wt") for cell_type in set(barcodes.values())}

    with gzip.open(forward_fq, "rt") as fq_f, gzip.open(reverse_fq, "rt") as fq_r:
        while True:
            # Read a chunk of records
            block_f = [fq_f.readline().strip() for _ in range(4 * chunk_size)]
            if not block_f[0]: # If the first element is an empty string (end of file)
                break

            block_r = [fq_r.readline().strip() for _ in range(4 * chunk_size)]

            for i in range(chunk_size):
                index = i * 4
                if not block_f[index]:
                    break  # Exit if EOF is reached within the chunk

                # Extract information for each record in the chunk
                bc = block_f[index].split()[0].split(":CR_")[1]
                cell_type = barcodes.get(bc)
                if cell_type:
                # Write to the appropriate output files based on barcode
                    forward_handles[cell_type].write(f"{block_f[index]}\n{block_f[index + 1]}\n+{block_f[index][1:]}\n{block_f[index + 3]}\n")
                    reverse_handles[cell_type].write(f"{block_r[index]}\n{block_r[index + 1]}\n+{block_r[index][1:]}\n{block_r[index + 3]}\n")

    # Close all file handles
    for handle in forward_handles.values():
        handle.close()
    for handle in reverse_handles.values():
        handle.close()

fastq_dir = sys.argv[1] # "/media/leon/Masha/ATAC/fastqs/SRR14048779"
sample = sys.argv[2] # "SRR14048779"
output_dir = sys.argv[3] # "/media/leon/Masha/ATAC/fastqs/SRR14048779/split"
barcodes_csv = sys.argv[4] # "/media/leon/Masha/ATAC/ATAC_barcodes/SRR14048779_barcodes.csv"

split_fastq(fastq_dir, sample, output_dir, barcodes_csv)