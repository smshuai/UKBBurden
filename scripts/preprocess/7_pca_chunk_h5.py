# conda activate cyvcf
import sys
import numpy as np
import h5py
from cyvcf2 import VCF

def median_impute(arr):
    ''' Impute missing GT (3) with median of the rest
    '''
    arr[arr==3] = np.int32(np.median(arr[arr<3]))
    return(arr)


def gt_counts(arr):
    unique, counts = np.unique(arr, return_counts=True) 
    print(unique)
    print(counts)

    
def read_vcf_into_arr(vcf_path):
    vcf_reader = VCF(vcf_path, gts012=True, threads=4)
    common_snvs = list()
    rare_snvs = list()
    for record in vcf_reader:
        if record.aaf >= 0.01 and record.aaf <= 0.99:
            # Common
            gt_arr = median_impute(record.gt_types.copy())
            common_snvs.append(gt_arr.copy())
        else:
            # Rare
            mac = (record.num_het + record.num_hom_alt * 2) if record.aaf < 0.01 else (record.num_het + record.num_hom_ref * 2)
            if mac >= 5:
                gt_arr = median_impute(record.gt_types.copy())
                rare_snvs.append(gt_arr)
    return(common_snvs, rare_snvs)


if __name__ == '__main__':
    vcf_list = np.loadtxt(sys.argv[1], str)
    # write with h5py
    with h5py.File(sys.argv[2], 'w') as f:
        first = True
        for vcf_path in vcf_list:
            print(vcf_path)
            common_snvs, rare_snvs = read_vcf_into_arr(vcf_path)
            print(len(common_snvs))
            print(len(rare_snvs))
            if first:
                first = False
                f.create_dataset("rare", data=np.array(rare_snvs), maxshape=(None, 200643))
                f.create_dataset("common", data=np.array(common_snvs), maxshape=(None, 200643))
            else:
                # append
                f['rare'].resize(f['rare'].shape[0] + len(rare_snvs), axis=0)
                f['rare'][-len(rare_snvs):] = np.array(rare_snvs)
                f['common'].resize(f['common'].shape[0] + len(common_snvs), axis=0)
                f['common'][-len(common_snvs):] = np.array(common_snvs)