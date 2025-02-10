import os
import sys
import scanpy as sc
import pandas as pd

adata_path = "/media/leon/Masha/ATAC/pipeline_test/processed_h5ads/snATAC_all_5kb.h5ad"
output_dir = "/media/leon/Masha/ATAC/ATAC_barcodes/"

cell_types = {
    '0': 'alpha', '1': 'beta', '2': 'beta', '3': 'beta', '4': 'alpha',
    '5': 'beta', '6': 'delta', '7': 'acinar', '8': 'alpha', '9': 'ductal',
    '10': 'gamma', '11': 'stellate', '12': 'immune', '13': 'EC'
}

ids_dct = {'JYH792': 'SRR14048778', 'MM110': 'SRR14048782', 'MM123': 'SRR14048753', 'MM124': 'SRR14048754', 'MM56': 'SRR14048758', 'MM59': 'SRR14048760', 'MM80': 'SRR14048766', 'MM86': 'SRR14048768', 'MM89': 'SRR14048771', 'MM95': 'SRR14048774', 'MM98': 'SRR14048777', 'MM108': 'SRR14048780', 'MM55': 'SRR14048757', 'MM61': 'SRR14048762', 'MM77': 'SRR14048763', 'MM78': 'SRR14048764', 'MM87': 'SRR14048769', 'MM93': 'SRR14048772', 'MM96': 'SRR14048775', 'JYH809': 'SRR14048779', 'MM109': 'SRR14048781', 'MM12': 'SRR14048783', 'MM120': 'SRR14048750', 'MM121': 'SRR14048751', 'MM122': 'SRR14048752', 'MM51': 'SRR14048755', 'MM54': 'SRR14048756', 'MM57': 'SRR14048759', 'MM60': 'SRR14048761', 'MM79': 'SRR14048765', 'MM81': 'SRR14048767', 'MM88': 'SRR14048770', 'MM94': 'SRR14048773', 'MM97': 'SRR14048776'}

def return_rev_compl(bc):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(bc))

adata = sc.read_h5ad(adata_path)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

samples = list(ids_dct.values())
for sample in samples:
    df = pd.DataFrame()
    df['barcode'] = adata.obs.index.str.split('_').str[1]
    if sample in ["SRR14048778", "SRR14048779", "SRR14048783"]:
        df['barcode'] = df['barcode'].map(return_rev_compl)
    df['cell_type'] = adata.obs.leiden.values.map(cell_types)
    df['donor'] = adata.obs.index.str.split('_').str[0]
    df['sample'] = df.donor.map(ids_dct)
    df = df[df['sample'] == sample]
    df = df[['barcode', 'cell_type']]
    name = f"{output_dir}/{sample}_barcodes.csv"
    df.to_csv(name, index=False)
