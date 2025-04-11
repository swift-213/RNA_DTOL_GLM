#!/usr/bin/env python

# script to count reads for each gene and each transcript

import argparse, sys, gzip

import pysam

from collections import defaultdict

import pandas as pd

import copy

def parseGff3(gff):
    #little function to parse the info line
    output = defaultdict(dict)
    for line in gff:
        if line == "" or line[0] == "#": continue
        
        gffObjects = line.strip().split("\t")
        
        chrom, feature, start, end, strand = gffObjects[0], gffObjects[2], int(gffObjects[3]), int(gffObjects[4]), gffObjects[6]
        
        try: info = dict([x.strip().split("=") for x in gffObjects[-1].strip(";").split(";")])
        except: raise ValueError("Problem parsing info section: " + gffObjects[-1])        
        
        if feature == "gene" or feature == "ncRNA_gene":
            output[chrom][info["ID"]] = {"type":feature, 'start':start, 'end':end, 'strand':strand, 'n_transcripts':0, 'n_aligned_transcripts':0, 'transcripts':{}}
        
        if feature == "mRNA":
            gene = info["Parent"]
            output[chrom][gene]['n_transcripts'] += 1
            output[chrom][gene]['transcripts'][info["ID"]] = {"type":feature, 'start':start, 'end':end, 'strand':strand, 'n_exons':0, 'n_CDSs':0,
                                                         'fivePrimeUTRs':[], 'exons':[], 'CDSs':[], 'threePrimeUTRs':[]}
        if feature in ("snRNA", "lnc_RNA", "ncRNA", "rRNA", "tRNA", "snoRNA", "scRNA", "transcript"):
            gene = info["Parent"]
            output[chrom][gene]['n_transcripts'] += 1
            output[chrom][gene]['transcripts'][info["ID"]] = {"type":feature, 'start':start, 'end':end, 'strand':strand, 'n_exons':0, 'exons':[]}
        
        if feature == "five_prime_UTR":
            transcript = info["Parent"]
            output[chrom][gene]['transcripts'][transcript]['fivePrimeUTRs'].append((start,end,))
        if feature == "three_prime_UTR":
            transcript = info["Parent"]
            output[chrom][gene]['transcripts'][transcript]['threePrimeUTRs'].append((start,end,))
        if feature == "exon":
            transcript = info["Parent"]
            output[chrom][gene]['transcripts'][transcript]['n_exons'] += 1
            output[chrom][gene]['transcripts'][transcript]['exons'].append((start,end,))
        if feature == "CDS":
            transcript = info["Parent"]
            output[chrom][gene]['transcripts'][transcript]['n_CDSs'] += 1
            output[chrom][gene]['transcripts'][transcript]['CDSs'].append((start,end,))
    return(output)


def get_count_data(geneData, bams, chrom_name_dict):
    chrom_names = bams[0].references
    for chrom in chrom_names:
        print(f"Doing chromosome {chrom}.", file=sys.stderr)
        if chrom_name_dict:
            try: chrom = chrom_name_dict[chrom]
            except:
                print(f"Warning: chromosome '{chrom}' not included in chromosome conversion file. Skipping.", file=sys.stderr)
                continue
        print(f"{len(geneData[chrom])} genes found.", file=sys.stderr)
        for gene_num, gene in enumerate(geneData[chrom].keys()):
            info = geneData[chrom][gene]
            counts = [bam.count(contig=chrom, start=info["start"]-1, stop=info["end"]) for bam in bams]
            yield [chrom, info["start"], info["end"], info["type"], info["strand"], gene, "NA"] + counts
            
            for transcript_num, transcript in enumerate(geneData[chrom][gene]["transcripts"].keys()):
                info = geneData[chrom][gene]["transcripts"][transcript]
                transcript_counts = [bam.count(contig=chrom, start=info["start"]-1, stop=info["end"]) for bam in bams]
                yield [chrom, info["start"], info["end"], info["type"], info["strand"], transcript, gene] + transcript_counts
                exon_info = geneData[chrom][gene]["transcripts"][transcript]['exons']
                for exon_start, exon_end in exon_info:
                    exon_counts = [bam.count(contig=chrom, start=exon_start-1, stop=exon_end) for bam in bams]
                    yield [chrom, exon_start, exon_end, 'exon', info["strand"], transcript, gene] + exon_counts

#i need to edit this code so that it counts the reads that are aligned to each aligned exon in the gene and then takes this count as the per gene count
#geneData[chrom][gene]['transcripts']

def subset_gff(df_of_aligned_exons, gff_object, alt_or_ref):
    prefix_gene_id = alt_or_ref + '_gene_id'
    if alt_or_ref == 'alt':
        chromosome_prefix = 'alt_chrom'
    else:
        chromosome_prefix = 'ref_chrom_num'
    reference_genes_with_exons = {prefix_gene_id : exons for prefix_gene_id, exons in df_of_aligned_exons.groupby(prefix_gene_id)}
    reference_chromosomes_with_genes = {chromosome : exons for chromosome, exons in df_of_aligned_exons.groupby(chromosome_prefix)}
    #get a dict where you have the gene_ID with each of its ref gene id's 
    refchrom_geneID = {chromosome: set(exons[prefix_gene_id])for chromosome, exons in reference_chromosomes_with_genes.items()}
    #we then create a new dict object that we will subset the old gff object into
    #we then subset so we only have reference genes that have aligned exons
    filter_by_gene = defaultdict(dict)
    for chromosome, aligned_gene_IDs in refchrom_geneID.items():
        for aliged_id in aligned_gene_IDs:
            filter_by_gene[str(chromosome)][aliged_id] = copy.deepcopy(gff_object[str(chromosome)][aliged_id])
    filter_by_transcript = defaultdict(dict)
    for ref_gene_id, exons in reference_genes_with_exons.items():
        transcript_IDs = set(exons[(alt_or_ref + '_transcript_id')])
        filter_by_transcript[str(exons[chromosome_prefix].iloc[0])][ref_gene_id] = copy.deepcopy(filter_by_gene[str(exons[chromosome_prefix].iloc[0])][ref_gene_id])
        filter_by_transcript[str(exons[chromosome_prefix].iloc[0])][ref_gene_id]['transcripts'] = {}
        filter_by_transcript[str(exons[chromosome_prefix].iloc[0])][ref_gene_id]['n_aligned_transcripts'] += len(transcript_IDs)
        #if wanted to only keep genes where all exons align can easily add in line of code here 
        for id in transcript_IDs:
            filter_by_transcript[str(exons[chromosome_prefix].iloc[0])][ref_gene_id]['transcripts'][id] = filter_by_gene[str(exons[chromosome_prefix].iloc[0])][ref_gene_id]['transcripts'][id]
    #we then take the new only aligned gene gff object and we go through and only keep exons that are perfectly aligned between the two 
    for reference_gene_id, exons in reference_genes_with_exons.items():
        transcript_IDs = set(exons[(alt_or_ref + '_transcript_id')])
        for transcript in transcript_IDs:
            exons_in_gff = filter_by_transcript[str(exons[chromosome_prefix].iloc[0])][reference_gene_id]['transcripts'][transcript]['exons']
            checked_exons=[]
            for start, end in exons_in_gff:
                if any(start == exons[(alt_or_ref +'_start')]) and any(end == exons[alt_or_ref + '_end']):
                    checked_exons.append((start, end))
            if checked_exons == exons_in_gff:
                continue
            else:
                filter_by_transcript[str(exons[chromosome_prefix].iloc[0])][exons[prefix_gene_id].iloc[0]]['transcripts'][transcript]['exons'] = checked_exons
    return filter_by_transcript





parser=argparse.ArgumentParser()

parser.add_argument("-b", "--bams", help="Input bam file(s) separated by spaces", action = "store", nargs = "+", required = True)
parser.add_argument("-g", "--gff", help="Gff3 annotation file file", action='store', required = True)
parser.add_argument("-ae", "--aligned_exon_csv", help="CSV file of aligned exon positions", action='store', required = True)
parser.add_argument("-r_o_a", "--ref_or_alt_assembly", help="Are we using the reference or alternate gff?", action='store', required = True)
parser.add_argument("-o", "--out_file", help="Output file name", action = "store")
parser.add_argument("--bam_to_gff_chrom_file", help="File relating chrom name in bam (1st column) to gff (2nd column)", action = "store")
parser.add_argument("--gff_to_bam_chrom_file", help="File relating chrom name in gff (1st column) to bam (2nd column)", action = "store")

args = parser.parse_args()

aligned_exons_df = pd.read_csv(args.aligned_exon_csv)
#aligned_exons_df = pd.read_csv('/home/s1929681/Desktop/filtered_aligned_exons.csv')

#"/data/martin/genomics/analyses/DTOL_insect_indels/whole_genome_files/VCFs/iyBomPrat1.1.vcf.gz"
print("Loading gene data...", file=sys.stderr)
with gzip.open(args.gff,"rt") if args.gff.endswith(".gz") else open(args.gff, "rt") as gff:
    geneData_unfiltered = parseGff3(gff)
    if args.ref_or_alt_assembly == 'alt':
        print('Using the alternate gff')
        geneData = subset_gff(aligned_exons_df, geneData_unfiltered, 'alt')
    elif args.ref_or_alt_assembly == 'ref':
        print('Using the reference gff')
        geneData = subset_gff(aligned_exons_df, geneData_unfiltered, 'ref') 

if args.gff_to_bam_chrom_file:
    with open(args.gff_to_bam_chrom_file, "rt") as chromfile:
        chrom_name_dict = dict([(line.split()[:2]) for line in chromfile])

if args.bam_to_gff_chrom_file:
    with open(args.bam_to_gff_chrom_file, "rt") as chromfile:
        chrom_name_dict = dict([(line.split()[2::-1]) for line in chromfile])
else:
    chrom_name_dict = None

#with gzip.open('/media/s1929681/Seagate/Frankie_DTOL_lep_project/phase_2_gff/ref_gff/iyBomPrat1.1_GCA_930367275.1.gff3',"rt") if '/media/s1929681/Seagate/Frankie_DTOL_lep_project/phase_2_gff/ref_gff/iyBomPrat1.1_GCA_930367275.1.gff3'.endswith(".gz") else open('/media/s1929681/Seagate/Frankie_DTOL_lep_project/phase_2_gff/ref_gff/iyBomPrat1.1_GCA_930367275.1.gff3', "rt") as gff:
#    geneData_unfiltered = parseGff3(gff)
#    geneData = subset_gff(aligned_exons_df, geneData_unfiltered, 'ref')

#with gzip.open("/media/s1929681/Seagate/Frankie_DTOL_lep_project/phase_2_gff/alt_gff/iyBomPrat1.1_alt.gff3","rt") if '/media/s1929681/Seagate/Frankie_DTOL_lep_project/phase_2_gff/alt_gff/iyBomPrat1.1_alt.gff3'.endswith(".gz") else open('/media/s1929681/Seagate/Frankie_DTOL_lep_project/phase_2_gff/alt_gff/iyBomPrat1.1_alt.gff3', "rt") as gff:
#    geneData = parseGff3(gff)



#bam reader
bams = [pysam.AlignmentFile(b, "rb") for b in args.bams]
#bams = pysam.AlignmentFile("/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/phase_bam_iyBomPrat1.1_GCA_930367275.1_1.bam", "rb")

if args.out_file:
    out = gzip.open(args.out_file,"wt") if args.out_file.endswith(".gz") else open(args.out_file, "wt")
else:
    out = sys.stdout

out.write("Scaffold\tStart\tEnd\tFeature\tStrand\tID\tParent\t")
out.write("\t".join(args.bams) + "\n")

print("Getting count data...", file=sys.stderr)
for data in get_count_data(geneData, bams, chrom_name_dict):
    out.write("\t".join([str(x) for x in data]) + "\n")

for bam in bams: bam.close()
out.close()

