import os
import sys
import gzip
import pandas as pd

def read_barcodes(file_path):
    df = pd.read_csv(file_path)
    barcodes = df.groupby('cell_type')['barcode'].apply(lambda x: set(x)).to_dict()
    return barcodes

def split_fastq(fastq_dir, sample, output_dir, barcodes_csv, chunk_size=100000):
    barcodes = read_barcodes(barcodes_csv)
    forward_fq = f"{fastq_dir}/{sample}/{sample}_S1_L001_R1_001.fastq.gz"
    reverse_fq = f"{fastq_dir}/{sample}/{sample}_S1_L001_R3_001.fastq.gz"
    barcode_fq = f"{fastq_dir}/{sample}/{sample}_S1_L001_R2_001.fastq.gz"
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    forward_handles = {cell_type: gzip.open(f"{output_dir}/{sample}_{cell_type}_1.fastq.gz", "wt") for cell_type in barcodes}
    reverse_handles = {cell_type: gzip.open(f"{output_dir}/{sample}_{cell_type}_2.fastq.gz", "wt") for cell_type in barcodes}

    with gzip.open(forward_fq, "rt") as fq_f, gzip.open(reverse_fq, "rt") as fq_r, gzip.open(barcode_fq, "rt") as fq_b:
        while True:
            # Read a chunk of records
            block_b = [fq_b.readline().strip() for _ in range(4 * chunk_size)]
            if not block_b[0]: # If the first element is an empty string (end of file)
                break

            block_f = [fq_f.readline().strip() for _ in range(4 * chunk_size)]
            block_r = [fq_r.readline().strip() for _ in range(4 * chunk_size)]

            for i in range(chunk_size):
                index = i * 4
                if not block_b[index]:
                    break  # Exit if EOF is reached within the chunk

                # Extract information for each record in the chunk
                seq_b = block_b[index + 1]
                new_header_f = f"{block_f[index].split()[0]}_{seq_b} {block_f[index].split()[1]}"
                new_header_r = f"{block_r[index].split()[0]}_{seq_b} {block_r[index].split()[1]}"

                # Write to the appropriate output files based on barcode
                for cell_type, barcode_list in barcodes.items():
                    if seq_b in barcode_list:
                        forward_handles[cell_type].write(f"{new_header_f}\n{block_f[index + 1]}\n+\n{block_f[index + 3]}\n")
                        reverse_handles[cell_type].write(f"{new_header_r}\n{block_r[index + 1]}\n+\n{block_r[index + 3]}\n")

    # Close all file handles
    for handle in forward_handles.values():
        handle.close()
    for handle in reverse_handles.values():
        handle.close()

fastq_dir = sys.argv[1] # "/media/leon/Masha/ATAC/fastqs"
sample = sys.argv[2] # [f"SRR140487{i + 50}" for i in range(34)]
output_dir = sys.argv[3] # "/media/leon/Masha/ATAC/cell_type_specific_fastqs"
barcodes_csv = sys.argv[4] # [f"/media/leon/Masha/ATAC/ATAC_barcodes/{sample}_barcodes.csv" for sample in SRR_ids]

split_fastq(fastq_dir, sample, output_dir, barcodes_csv)