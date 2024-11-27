import math
import numpy as np
import argparse
import subprocess
from cyvcf2 import Writer, VCF
from multiprocessing import Pool

def filter_vcf(out, vcf_path, snp = 0, biallelic = 0, indel = 0, DP = 0, GQ = 0, AD = 0,
               covg_file = "", frac_called = 0, quant=0.9, region = '', samples=''):
    """
    Default no filters applied.
    Needs an "average_depth.txt" file in same dir which specifies average depth of each isolate in the vcf.
    """

    print(f"Initiating for {region}.")
    vcf = VCF(vcf_path)
    
    sample_indices = []
    if samples:
        assert isinstance(samples, list), "Samples input must be a list"
        sample_indices = [i for i, j in enumerate(vcf.samples) if j in samples]
        vcf.set_samples(samples)
    num_samples = len(vcf.samples)
    quant_x = math.ceil(quant * num_samples)

    # Initializing out-VCF
    w = Writer(out, vcf)
    
    # Load average read depth for isolates
    np_avgs = np.zeros(num_samples)
    if covg_file:
        avgs_array = []
        with open(covg_file,"r") as file:
            for line in file:
                avgs_array.append(float(line.strip()))
        np_avgs = np.array(avgs_array)
    np_max = np_avgs + 4 * np.sqrt(np_avgs) # need to define whether or not coverage file is provided
    if sample_indices != []:
        np_max = np_max[np.array(sample_indices)]
    
    # Filters
    def check_snp(row):
        return row.is_snp
    def check_biallelic(row):
        return len(row.ALT) == 1
    def num_that_pass_metrics(row):
        passed = 0
        counter = -1
        gqs = ['GQ', 'RGQ']['RGQ' in row.FORMAT]
        gqs = row.format(gqs).flatten()
        dps = row.format('DP').flatten()
        for gq, dp in zip(gqs, dps):
            counter += 1
            if gq >= GQ and dp >= DP and row.gt_depths[counter] < np_max[counter]:
                passed += 1
        return passed >= quant_x
    def check_called_frac(row):
        return row.num_called >= frac_called*len(vcf.samples)
    def check_indel(row):
        return row.is_indel
    def check_not_indel(row):
        return not row.is_indel
    def check_mt_minus(row):
        start, end = 465219, 679919
        return row.POS < start or row.POS > end

    # Assemble filters
    condition_checks = []
    if snp:
        condition_checks.append(check_snp)
    if biallelic:
        condition_checks.append(check_biallelic)
    if indel == 1:
        condition_checks.append(check_indel)
    elif indel == 2:
        condition_checks.append(check_not_indel)
    # if DP or GQ or AD or covg_file:
    if GQ or DP or AD or covg_file:
        condition_checks.append(num_that_pass_metrics)
    if frac_called:
        condition_checks.append(check_called_frac)
    if region == 'chromosome_06':
        condition_checks.append(check_mt_minus)
    if len(condition_checks) == 0:
        raise Exception("Please specify at least one filtering condition.")

    # Apply filters
    for row in vcf(region):
        if all(condition(row) for condition in condition_checks):
            w.write_record(row)
    print(f"Done for {region}. Closing this file.")
    w.close(), vcf.close()

def make_header(VCF, out_name):
    with open(out_name, "w") as file:
        file.write(VCF.raw_header)

def parallelize_cyvcf2(arg_list: tuple, chromosomes: list, processes: int, out_name: str):
    """
    Parallelize cyvcf2. Specify processes as '-1' to allocate one process per chromosome.
    """
    
    if processes == -1:
        processes = len(chromosomes)
        
    in_vcf = VCF(arg_list[0])
    if arg_list[-1] != []:
        in_vcf.set_samples(arg_list[-1])
    
    make_header(in_vcf, out_name)
    
    temp_name = 'temp_dir_0606'
    try:
        subprocess.run(['mkdir', temp_name])
    except:
        subprocess.run(['rm', '-rf', temp_name])
        subprocess.run(['mkdir', temp_name])

    inputs = [(f'{temp_name}/{i}.vcf',)     # chromosome-specific VCF name
              + arg_list[:-2]               # arguments until 'region' and 'samples'
              + (i,)                        # region
              + (arg_list[-1],)             # samples
              for i in chromosomes]        

    # Run parallelized computation
    if __name__ == '__main__':
        with Pool(processes) as p:
            p.starmap(filter_vcf, inputs)
    
    # Sorting chromosome-specific VCF filenames
    chromosome_files = sorted([f"{temp_name}/chromosome_{i:02d}.vcf" for i in range(1, 18)])
    
    # Merge sorted VCFs into final VCF
    with open(out_name, 'a') as a_file:
        for chromosome_file in chromosome_files:
            with open(chromosome_file, "r") as chr_file:
                line = chr_file.readline()
                while line.startswith('#'):
                    line = chr_file.readline()  # Skip headers after first file
                a_file.write(line)
                for line in chr_file:
                    a_file.write(line)
                    
    subprocess.run(['rm', '-rf', temp_name])

def main():
    # Parsing input arguments
    parser = argparse.ArgumentParser(description="Filter a VCF file using various conditions.")
    parser.add_argument("--out", required=True, help="Output VCF file.")
    parser.add_argument("--vcf_path", required=True, help="Input VCF file path.")
    parser.add_argument("--snp", type=int, default=0, help="Filter SNPs (1 to enable).")
    parser.add_argument("--biallelic", type=int, default=0, help="Filter biallelic sites (1 to enable).")
    parser.add_argument("--indel", type=int, default=0, help="Filter indels (1: include, 2: exclude).")
    parser.add_argument("--DP", type=int, default=0, help="Minimum depth.")
    parser.add_argument("--GQ", type=int, default=0, help="Minimum genotype quality.")
    parser.add_argument("--AD", type=int, default=0, help="Minimum allele depth.")
    parser.add_argument("--covg_file", default="", help="File with average depth per isolate.")
    parser.add_argument("--frac_called", type=float, default=0, help="Fraction of samples with calls.")
    parser.add_argument("--quant", type=float, default=0.9, help="Proportion of samples meeting criteria.")
    parser.add_argument("--chromosomes", default='', help="Regions to filter on from file (e.g., chromosome_01:10000-100000).")
    parser.add_argument("--samples", nargs='*', default=[], help="List of sample IDs.")
    parser.add_argument("--processes", type=int, default=-1, help="Number of parallelizations to run.")
    args = parser.parse_args()

    # Printing to terminal the options chosen.
    if args.snp:
        print("Must be a SNP")
    if args.biallelic:
        print("Must be biallelic")
    if args.indel == 1:
        print("Must be an indel")
    elif args.indel == 2:
        print("Must not be an indel")
    if args.GQ or args.DP or args.AD or args.covg_file:
        print(f"Must have a DP of {str(args.DP)}, a GQ of {str(args.GQ)}, and an AD of {str(args.AD)}.")
    if args.covg_file:
        print("Must do the coverage thing.")
    if args.frac_called:
        print("Fraction of sites called must be " + str(args.frac_called))
    if args.samples:
        print(f"Running on the following samples: {args.samples}")
    statement = [f"Running {args.processes} parallelizations.", "Running one parallel process per chromosome."][args.processes==-1]
    print(statement)

    # Executing parallelize_cyvcf2
    chromosomes = []
    with open(args.chromosomes, "r") as file:
        for line in file:
            chromosomes.append(line.strip())
    c =       (args.vcf_path, 
               args.snp, 
               args.biallelic, 
               args.indel,
               args.DP, args.GQ, args.AD, 
               args.covg_file, 
               args.frac_called,
               args.quant, 
               args.chromosomes, 
               args.samples,
              )
    parallelize_cyvcf2(c, chromosomes, args.processes, args.out)

    print("Finished.")

if __name__ == "__main__":
    main()
