#!/usr/bin/env python

# script to split a bam file into two based on phased haplotypes.
# Only reads that perfectly match one or the other haplotype will be kept
# Option to use REF and ALT bases instead for caes of an alt assembly aligned to ref. In these cases 

import argparse, sys, gzip, os

from collections import defaultdict

import numpy as np

import pysam

import cyvcf2

import tqdm as tqdm

#Gets all the SNPs in the windows?
def get_phased_variants(vcf, chrom, start, end, sample_idx=None, useREFandALT=False):
    variants = {}
    for variant in vcf(chrom+":"+str(start)+"-"+str(end)):
        variant_alleles = [variant.REF] + variant.ALT
        if len(variant_alleles) < 2: continue
        if not useREFandALT:
            #set alleles based on an individual's phased genotype
            geno = variant.genotypes[sample_idx if sample_idx else 0][:-1]
            if len(geno) < 2 or len(set(geno))==1: continue
            variant_alleles = [variant_alleles[i] for i in geno]
        
        variants[variant.POS] = variant_alleles
    
    return variants


def get_phase_for_reads(reads, variants, min_matches=1, min_base_qual=0, sample_idx=None, useREFandALT=False, chrom_name_dict=None):
    
    #give up if not enough informative variants
    if len(variants) < min_matches: return None 
    
    #number of matches to haplotypes 0 and 1, summed across all reads
    matches = [0,0]
    for read in reads:
        for read_pos, ref_pos in read.get_aligned_pairs():
            if read_pos != None and ref_pos != None and ref_pos+1 in variants and qualDict[read.qual[read_pos]] >= min_base_qual:
                #does it match the first or second haplotype (or neither)?
                try:
                    hap_match = variants[ref_pos+1].index(read.query_sequence[read_pos])
                    matches[hap_match] += 1
                except:
                    continue
    
    if min(matches) == 0 and max(matches) >= min_matches:
        return 0 if matches[0] > 0 else 1
    else:
        return None


qualDict = dict(zip(("!", "\"", "#", "$", "%", "&",
                     "'", "(", ")", "*", "+", ",",
                     "-", ".", "/", "0", "1", "2",
                     "3", "4", "5", "6", "7", "8",
                     "9", ":", ";", "<", "=", ">",
                     "?", "@", "A", "B", "C", "D",
                     "E", "F", "G", "H", "I", "J", "K"), range(43)))


parser=argparse.ArgumentParser()

parser.add_argument("-b", "--bam", help="Input bam file", action = "store", required = True)
parser.add_argument("-v", "--vcf", help="Completely phased VCF file", action='store', required = True)
parser.add_argument("-o", "--out_prefix", help="Name prefix for output bams", action = "store", required = True)
parser.add_argument("--tmpdir", help="Temporary directory for bam file sorting", action = "store", default=".")
parser.add_argument("--samtools_threads", help="Threads for samtools sorting steps", action = "store", type=int, default=1)
parser.add_argument("--use_REF_and_ALT", help="Use REF allele as first haplotype and ALT allele as second", action = "store_true")
parser.add_argument("--sampleID", help="Sample ID if not using first sample in VCF", action = "store")
parser.add_argument("--min_map_qual", help="Minumum mapping quality for reads", action = "store", type=int, default=0)
parser.add_argument("--min_base_qual", help="Minumum base quality for matching haplotypes", action = "store", type=int, default=20)
parser.add_argument("--min_matches", help="Minumum high quality bases to match diagnostic allele (summed across read pair)", action = "store", type=int, default=1)
parser.add_argument("--bam_to_vcf_chrom_file", help="File relating chrom name in bam (1st column) to vcf (2nd column)", action = "store")
parser.add_argument("--vcf_to_bam_chrom_file", help="File relating chrom name in vcf (1st column) to bam (2nd column)", action = "store")
parser.add_argument("--max_pair_dist", help="Maximum distance to look for paired reads (must be less than window size)", action = "store", type=int, default=100000)
parser.add_argument("--run_quietly", help="Prevents printing of every window", action = "store", default=False, required = False)


#args = parser.parse_args("-b temp.bam -v /data/martin/genomics/analyses/DTOL_insect_indels/whole_genome_files/VCFs/iyBomPrat1.1.vcf.gz --use_REF_and_ALT -o test --vcf_to_bam_chrom_file iyBomPrat1.1_chromosomes.txt".split())

args = parser.parse_args()
#vcf = cyvcf2.VCF('/Volumes/Seagate/Frankie_DTOL_lep_project/outputs/VCF/ilEupSimi1.1.vcf.gz')
#vcf = cyvcf2.VCF('/Volumes/Seagate/Frankie_DTOL_lep_project/outputs/VCF/iyBomPrat1.1.vcf.gz')
#"/data/martin/genomics/analyses/DTOL_insect_indels/whole_genome_files/VCFs/iyBomPrat1.1.vcf.gz"
vcf = cyvcf2.VCF(args.vcf, samples= [args.sampleID] if args.sampleID else None)

sample_idx = None if args.use_REF_and_ALT else 0

#with open('/Users/frankieswift/OneDrive/RA_Work/Indel_Project/data_set_chroms/ilEupSimi1.1_chromosomes.txt', "rt") as chromfile:
#        chrom_name_dict = dict([(line.split()[2::-1]) for line in chromfile])


#with open('/Users/frankieswift/OneDrive/RA_Work/Indel_Project/data_set_chroms/iyBomPrat1.1_chromosomes.txt', "rt") as chromfile:
#        chrom_name_dict = dict([(line.split()[2::-1]) for line in chromfile])

if args.bam_to_vcf_chrom_file:
    with open(args.bam_to_vcf_chrom_file, "rt") as chromfile:
        chrom_name_dict = dict([(line.split()[:2]) for line in chromfile])

if args.vcf_to_bam_chrom_file:
    with open(args.vcf_to_bam_chrom_file, "rt") as chromfile:
        chrom_name_dict = dict([(line.split()[2::-1]) for line in chromfile])
else:
    chrom_name_dict = None


#bam = pysam.AlignmentFile("/data/martin/genomics/analyses/DTOL_insect_indels/RNA_seq/ERR7113577_Aligned.sortedByCoord.out.bam", "rb")
#bam = pysam.AlignmentFile("/Volumes/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/ERR7113577_Aligned.sortedByCoord.out.bam", "rb")
#bam = pysam.AlignmentFile("/Volumes/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/ERR6286704_Aligned.sortedByCoord.out.bam", "rb")
#bam reader
bam = pysam.AlignmentFile(args.bam, "rb")

chrom_names = bam.references
chrom_lengths = dict(zip(chrom_names, bam.lengths))

#output bams
outBams = [pysam.AlignmentFile(args.out_prefix + "_1_unsorted.bam", "wb", template = bam),
           pysam.AlignmentFile(args.out_prefix + "_2_unsorted.bam", "wb", template = bam)]


step_size = args.max_pair_dist

window_size = step_size*2

for chrom in chrom_names:
 
    if chrom_name_dict:
            try: vcf_chrom = chrom_name_dict[chrom]
            except:
                print(f"Warning: chromosome '{chrom}' not included in chromosome conversion file. Skipping.", file=sys.stderr)
                continue
    print(f"Analysing chromosome {chrom}", file=sys.stderr)
    previous_read_names = set()
    for window_end in range(window_size, chrom_lengths[chrom] + step_size, step_size):
        window_start = window_end - window_size
        
        if not args.run_quietly:
            print(chrom, window_start+1, window_end, file=sys.stderr)
        
        midpoint = window_start + step_size
        
        window_reads = list(bam.fetch(chrom, window_start, window_end))
        
        #mapping qual fiter
        window_reads = [read for read in window_reads if read.mapping_quality >= args.min_map_qual]
        
        #those in the first half of the window that were processed previously can be removed
        window_reads = [read for read in window_reads if read.reference_start > midpoint or read.query_name not in previous_read_names]
        
        #How many times the reads are present??
        pair_names, pair_counts = np.unique([read.query_name for read in window_reads], return_counts=True)
        
        name_counts = dict(zip(pair_names, pair_counts))
        
        reads_by_pair = defaultdict(list)
        
        for i, read in enumerate(window_reads):
            if read.reference_end < midpoint or name_counts[read.query_name] > 1:
                reads_by_pair[read.query_name].append(read)
        
        npairs = len(reads_by_pair)
        
        if not args.run_quietly:
            print(f"{npairs} read pairs.", file=sys.stderr)
        
        if npairs ==0:
            print("\n", file=sys.stderr)
            continue
        
        #get the start and end for variants, which will possibly extend a little either way because of the reads that overlap
        variants_start = min([read.reference_start for read in window_reads]) + 1 #also make it 1-based for vcf
        variants_end = max([read.reference_end for read in window_reads]) + 1
        
        variants = get_phased_variants(vcf, chrom = chrom_name_dict[chrom] if chrom_name_dict else chrom,
                                       start=variants_start, end=variants_end,
                                       sample_idx=sample_idx, useREFandALT=args.use_REF_and_ALT)

        if not args.run_quietly:
            print(f"{len(variants)} phased variants.", file=sys.stderr)

        pairs_phased = 0
        
        for reads in reads_by_pair.values():
            #phase = get_phase_for_reads(reads, variants, min_matches=1, min_base_qual=20)
            phase = get_phase_for_reads(reads, variants, min_matches=args.min_matches, min_base_qual=args.min_base_qual)
            
            if phase is not None:
                pairs_phased += 1
                for read in reads:
                    outBams[phase].write(read)
        
        if not args.run_quietly:
            print(f"{pairs_phased} read pairs were phased.\n", file=sys.stderr)
        
        previous_read_names = set(reads_by_pair.keys())



bam.close()
outBams[0].close()
outBams[1].close()


#sort the new bams by position and remove name_sorted_ones
pysam.sort("-@", str(args.samtools_threads), "-T", args.tmpdir, "-o", args.out_prefix + "_1.bam", args.out_prefix + "_1_unsorted.bam")
pysam.index(args.out_prefix + "_1.bam")
os.remove(args.out_prefix + "_1_unsorted.bam")

pysam.sort("-@", str(args.samtools_threads), "-T", args.tmpdir, "-o", args.out_prefix + "_2.bam", args.out_prefix + "_2_unsorted.bam")
pysam.index(args.out_prefix + "_2.bam")
os.remove(args.out_prefix + "_2_unsorted.bam")


sys.stderr.write("\nDone.\n")
