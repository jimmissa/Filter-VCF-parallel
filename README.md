This folder contains a command-line Pythonic program to filter a VCF file, using cyvcf2, while parallelizing filtration across each chromosome to dramatically increase the speed of filtration.

The program is contained in filter_vcf.py and can be run from the command-line as follows:

python filter_vcf.py --out filtered.vcf --vcf-path path/to/file.vcf --snp 1 --biallelic 1 --indel 2 --DP 3 --GQ 30 --AD 5 --covg_file average_depth_new.txt --frac_called 0.5 --quant 0.9 --chromosomes chromosomes.txt --samples ISOLATE1 ISOLATE2 ISOLATE3 --processes -1

A second way to run it is also made available via the run_filter_vcf.sh file. This file already sets up all the values for the input variables you may want, and then runs the Python program. You just need to call it from the command-line after setting it up: 

bash run_filter_vcf.sh 

All you need to do is edit the values of the input variables. If you are running from the command-line, you can exclude any input argument you do not want. If you are running directly running the run_filter_vcf.sh file, then just leave the filters you don't want to use as null values (0s for numerical filters like GQ, or an empty string if you don't want to specify a coverage file).

There are two additional files that may go alongside this program (other than the input VCF). One is a file that specifies the regions you want to parallelize filtration across. In this directory, it is called chromosomes.txt. It is newline-separated, specifying each region that filtration must be done on per line. This file is required. The second file is optional but allows for filtration to be done on in relation to the average sequencing depth of each individual in the VCF. For example, it may be desired to filter sites where many individuals have a sequencing depth above twice the average to exclude regions with copy number variation.
