#okay so the aim is to take the exon positions in both the ref and the alt and determine if they are actually aligned in the haplotypes and if they are aligned if the region is classes as a exon in the alt genome gff!
#So we will need the ref and alt bamfiles to ensure they are aligned as well as the ref and alt gffs!

#need to add in log file capabilities before cluster

import pysam
import gffpandas.gffpandas as gffpd
import pandas as pd
from tqdm import tqdm
import numpy as np
from argparse import ArgumentParser
from general.genfun import chromosome_dict_maker

#functions
def exon_gene_mrna_position_finder(input_df):
    exon_id=[]
    gene_ID = None  
    mRNA_ID = None
    for row, values in tqdm(input_df.iterrows()):
        if values['type']=='exon':
            if values['Parent'] == mRNA_ID:
                exon_id.append([values['seq_id'], values['type'], values['start'], values['end'], gene_ID, mRNA_ID, values['strand']])
                continue
            elif mRNA_ID == None:
                exon_id.append([values['seq_id'], values['type'], values['start'], values['end'], gene_ID, mRNA_ID, values['strand']])
                continue
        elif values['type'] == 'mRNA':
            if gene_ID == values['Parent']:
                mRNA_ID = values['ID']
                continue
        if values['type'] == 'gene':
            gene_ID = values['ID']
            continue
    return exon_id

def long_short_switch(chromosome):
    #print(chromosome)
    chromosome = str(chromosome)
    if chromosome in long_to_short:
        return long_to_short[chromosome]
    elif chromosome in short_to_long:
        return short_to_long[chromosome]
    return 'Na'


def gene_transcript_exon_info(gene_id_dict): 
    gene_info_dict = {}
    for gene, info in tqdm(gene_id_dict.items(), total = len(gene_id_dict)):
        gene_info = info.nunique()
        no_transcripts = int(gene_info['mRNA_ID'])
        if gene_info['start'] == gene_info['end']:
            number_exons = int(gene_info['start'])
        else:
            number_exons = len(set([tuple([row.start, row.end]) for row in info.itertuples(index = False)]))
        mRNA_ids = set([row.mRNA_ID for row in info.itertuples(index = False)])
        exon_positions = set([tuple([row.start, row.end]) for row in info.itertuples(index = False)])
        chromosome = set([row.chromosome for row in info.itertuples(index = False)])
        long_form_chrom = set([row.numeric_chrom for row in info.itertuples(index = False)])
        meta_info = [no_transcripts, number_exons, mRNA_ids,exon_positions, chromosome, long_form_chrom]
        gene_info_dict[gene] = meta_info
    return gene_info_dict

def appending_the_exon_df(all_aligned_alts, exisiting_list_name):
    for row in all_aligned_alts.itertuples(index = False):
        exisiting_list_name.append([values.numeric_chrom, values.chromosome, values.feature_type,
        values.start, values.end, values.gene_ID, values.mRNA_ID, values.strand,
        row.chromosome, row.feature_type, int(row.start), int(row.end),
        row.gene_ID, row.mRNA_ID, read.is_forward , row.strand])

def read_quaility_checks(reads, ref_chromosome_name_dict, alt_chromosome_name_dict, mapping_qual):
    passes = True
    if reads.reference_name not in ref_chromosome_name_dict:
        passes = False
    elif reads.mapping_quality != int(mapping_qual):
        passes = False
    if len(reads.query_name.split('|')) != 1:
        if reads.query_name.split('|')[2] not in alt_chromosome_name_dict:
        #elif reads.query_name not in alt_chromosome_name_dict:
            passes = False
    else:
        if reads.query_name not in alt_chromosome_name_dict:
            passes = False
    return passes

def multiple_genes_alinging_to_one(exons, string_gene_id,thrown_out_ids,all_gene_ids_with_the_aligned_ids, No_duplicate_start_end_positions, gene_id):
    all_alt_gene_IDs=set()
    for i in range(0, len(exons), 1):
        if exons[string_gene_id].iloc[i] not in all_alt_gene_IDs:
            if exons[string_gene_id].iloc[i] in all_gene_ids_with_the_aligned_ids.keys():
                if exons[string_gene_id].iloc[i] not in thrown_out_ids:
                    thrown_out_ids.add(exons[string_gene_id].iloc[i])
                    if all_gene_ids_with_the_aligned_ids[exons[string_gene_id].iloc[i]] in No_duplicate_start_end_positions:
                        del No_duplicate_start_end_positions[all_gene_ids_with_the_aligned_ids[exons[string_gene_id].iloc[i]]]
                        continue
                    else:
                        continue
                else:
                    continue
            else:
                all_gene_ids_with_the_aligned_ids[exons[string_gene_id].iloc[i]] = gene_id
                thrown_out_ids.add(exons[string_gene_id].iloc[i])
        else:
            continue

def single_gene_alinging(exons, string_gene_id, all_gene_ids_with_the_aligned_ids, thrown_out_ids, No_duplicate_start_end_positions, gene_id, multiple_chromosomes_aligning):
    break_out_of_this_iteration = False
    if exons[string_gene_id].iloc[0] in all_gene_ids_with_the_aligned_ids.keys():
        if exons[string_gene_id].iloc[0] not in thrown_out_ids:
            thrown_out_ids.add(exons[string_gene_id].iloc[0])
            if all_gene_ids_with_the_aligned_ids[exons[string_gene_id].iloc[0]] in No_duplicate_start_end_positions:
                del No_duplicate_start_end_positions[all_gene_ids_with_the_aligned_ids[exons[string_gene_id].iloc[0]]]
                break_out_of_this_iteration = True
            else:
                break_out_of_this_iteration = True
        else:
            break_out_of_this_iteration = True
    else:
        all_gene_ids_with_the_aligned_ids[exons[string_gene_id].iloc[0]] = gene_id
        if multiple_chromosomes_aligning == True:
            break_out_of_this_iteration = True
    return break_out_of_this_iteration

def filter_longest_exons(df, start_col, end_col):
    """
    Filters out shorter exons that share the same start or end positions.
    Keeps only the longest exon per unique start/end.
    """
    # Remove duplicates by keeping the longest range
    df = df[df[end_col] == df.groupby(start_col)[end_col].transform('max')]
    df = df[df[start_col] == df.groupby(end_col)[start_col].transform('min')]
    return df

def dict_values_to_df(dictionary):
    values_as_list=[]
    for key, value in dictionary.items():
        for row in value.itertuples():
            values_as_list.append(row)
    return pd.DataFrame(values_as_list)

def quaility_check_exons(gene_exon_dict, reference_or_alternate):
    No_duplicate_start_end_positions={}
    all_gene_ids_with_the_aligned_ids = {}
    thrown_out_ids = set()
    if reference_or_alternate != False:
        print('ref and alt dataframes')
        string_gene_id = str(reference_or_alternate) + '_gene_id'
        chromosome = str(reference_or_alternate) + '_chrom'
    else:
        print('Single dataframe')
        string_gene_id= 'gene_ID'
        chromosome = 'chromosome'
    for gene_id, exons in tqdm(gene_exon_dict.items(), total = len(gene_exon_dict)):
        if int(exons.nunique()[chromosome]) == 1:
            if int(exons.nunique()[string_gene_id]) == 1:
                if single_gene_alinging(exons, string_gene_id, all_gene_ids_with_the_aligned_ids, thrown_out_ids, No_duplicate_start_end_positions, gene_id, False) == True:
                    continue
                else:
                    pass
            else:
                multiple_genes_alinging_to_one(exons, string_gene_id,thrown_out_ids,all_gene_ids_with_the_aligned_ids, No_duplicate_start_end_positions , gene_id)
                continue
        else:
            if int(exons.nunique()[string_gene_id]) == 1:
                if single_gene_alinging(exons, string_gene_id, all_gene_ids_with_the_aligned_ids, thrown_out_ids, No_duplicate_start_end_positions, gene_id, True) == True:
                    continue
                else:
                    pass
            else:
                multiple_genes_alinging_to_one(exons, string_gene_id,thrown_out_ids,all_gene_ids_with_the_aligned_ids, No_duplicate_start_end_positions, gene_id)
                continue
        non_nested_exons = []
        for index, row in exons.iterrows():
            if reference_or_alternate != False:
                if not any((exons.ref_start < row.ref_start) & (exons.ref_end > row.ref_end) & (exons.index != index)):
                    non_nested_exons.append(index)
                if not any((exons.alt_start < row.alt_start) & (exons.alt_end > row.alt_end) & (exons.index != index)):
                    non_nested_exons.append(index)
            else:
                if not any((exons.start < row.start) & (exons.end > row.end) & (exons.index != index)):
                    non_nested_exons.append(index)
        #exon_counter_dict = counter(non_nested_exons)
        #all_non_nested_exons = removing_non_duplicated_values(non_nested_exons, exon_counter_dict)
        exons = (exons.loc[list(set(non_nested_exons))]).sort_index()
        #here we are trying to remove exons that have the same start or end position and taking only the longest one
        if reference_or_alternate != False:
            if any(len(exons) != exons[['ref_start', 'ref_end', 'alt_start', 'alt_end']].nunique()):
                # Apply filtering to remove redundant shorter exons
                ref_exons = filter_longest_exons(exons, 'ref_start', 'ref_end')
                ref_exons = filter_longest_exons(ref_exons, 'alt_start', 'alt_end')
                alt_exons = filter_longest_exons(exons, 'alt_start', 'alt_end')
                alt_exons = filter_longest_exons(alt_exons, 'ref_start', 'ref_end')
                # Ensure both reduced DataFrames match
                if ref_exons.sort_index().equals(alt_exons.sort_index()):
                    No_duplicate_start_end_positions[gene_id] = alt_exons
                else:
                    print('interesting')
            else:
                No_duplicate_start_end_positions[gene_id] = exons
        else:
            if any(len(exons) != exons[['start', 'end']].nunique()):
                # Apply filtering to remove redundant shorter exons
                ref_exons = filter_longest_exons(exons, 'start', 'end')
                No_duplicate_start_end_positions[gene_id] = ref_exons
            else:
                No_duplicate_start_end_positions[gene_id] = exons
    return No_duplicate_start_end_positions


############################################################################################################################################################# 

parser = ArgumentParser()
#inputs
parser.add_argument("-chromosomes", "--chromosome_name_files", help="File containing the chromosome names for the reference genome", required=True)
parser.add_argument("-bam", "--bamfile", help="bamfile for alignment between the reference and alternate file", required=True)
parser.add_argument("-alt_fai", "--alt_fasta_index", help="index file for alternate fasta file", required=True)
parser.add_argument("-ref_annotation", "--reference_annotation_file", help="annotation file for the reference genome", required=True)
parser.add_argument("-alt_annotation", "--alternate_annotation_file", help="annotation file for the alternate genome", required=True)

#outputs
parser.add_argument("-lf", "--log_file", action='store', help="logfile to output code response", required=False)
parser.add_argument("-aligned_exons_output", "--aligned_exons_output", action='store', help="csv file output of aligned exons between the reference and alternate genomes", required=True)
parser.add_argument("-exon_alignment_stats", "--exon_alignment_stats", action='store', help="csv file output for stats of alignment_of_exons", required=True)



args = parser.parse_args()


print('Loading in files')
chromosomes = args.chromosome_name_files
bamfile = pysam.AlignmentFile(args.bamfile, "rb")
ref_annotation = gffpd.read_gff3(args.reference_annotation_file)
alt_annotation = gffpd.read_gff3(args.alternate_annotation_file)
alt_fasta_index = args.alt_fasta_index



chromosomes = ("/home/fswift/0_scripts/0_git/ac3/Chromosome_name_conversion_files/ilAglIoxx1.1_chromosomes.txt")
#Reading in ref bam
bamfile = pysam.AlignmentFile("/mnt/loki/martin/frankie/RNA_seq_glm/1_analysis_files/0_sorted_bams/ilAglIoxx1.1_GCA_905147045.1_ref_num.sort.bam", "rb")
#Reading in ref gff
ref_annotation = gffpd.read_gff3('/mnt/loki/martin/frankie/RNA_seq_glm/raw_files/gff_files/ilAglIoxx1.1_GCA_905147045.1.gff3')
#Reading in alt gff
alt_annotation = gffpd.read_gff3('/mnt/loki/martin/frankie/RNA_seq_glm/raw_files/gff_files/ilAglIoxx1.1_GCA_905147125.1_liftover.gff3_polished')


#chromosomes = ("/home/s1929681/One_drive_file_copies_25_01_15/RA_work/Indel_Project/data_set_chroms/ilAmpTrag2.1_chromosomes.txt")
#Reading in ref bam
#bamfile = pysam.AlignmentFile("/media/s1929681/Seagate_B/Frankie_DTOL_lep_project/outputs/samtools/ilAmpTrag2.1_alignment.sort.bam", "rb")
#Reading in ref gff
#ref_annotation = gffpd.read_gff3('/media/s1929681/Seagate_B/Frankie_DTOL_lep_project/phase_2_gff/ref_gff/ilAmpTrag2.1_GCA_905220435.1.gff3')
#Reading in alt gff
#alt_annotation = gffpd.read_gff3('/media/s1929681/Seagate_B/Frankie_DTOL_lep_project/phase_2_gff/alt_gff/ilAmpTrag2.1_alt.gff3')
#alt_annotation = gffpd.read_gff3('/media/s1929681/Seagate_B/Frankie_DTOL_lep_project/liftoff/outputs/ilAmpTrag2.1_GCA_905220425.1_liftover.gff3_polished')

##iyBomPrat1.1
#chromosomes = ("/home/s1929681/One_drive_file_copies_25_01_15/RA_work/Indel_Project/data_set_chroms/iyBomPrat1.1_chromosomes.txt")
#Reading in ref bam
#bamfile = pysam.AlignmentFile("/media/s1929681/Seagate_B/Frankie_DTOL_lep_project/outputs/samtools/iyBomPrat1.1_alignment.sort.bam", "rb")
#Reading in ref gff
#ref_annotation = gffpd.read_gff3('/media/s1929681/Seagate_B/Frankie_DTOL_lep_project/phase_2_gff/ref_gff/iyBomPrat1.1_GCA_930367275.1.gff3')
#Reading in alt gff
#alt_annotation = gffpd.read_gff3('/media/s1929681/Seagate_B/Frankie_DTOL_lep_project/phase_2_gff/alt_gff/iyBomPrat1.1_alt.gff3')
#alt_annotation = gffpd.read_gff3('/media/s1929681/Seagate_B/Frankie_DTOL_lep_project/liftoff/outputs/iyBomPrat1.1_GCA_930367225.1_liftover.gff3')
#Getting the start and end positions for each chromosome in alt in order to do the array approach for SNP and monomorphic sites 
#with open("/media/s1929681/Seagate_B/Frankie_DTOL_lep_project/phase_two_fasta/ilAmpTrag2.1_GCA_905220425.1.fa.fai", "rt") as lf:
#with open("/media/s1929681/Seagate_B/Frankie_DTOL_lep_project/phase_two_fasta/iyBomPrat1.1_GCA_930367225.1.fa.fai", "rt") as lf:
with open("/mnt/loki/martin/frankie/RNA_seq_glm/raw_files/fasta_files/ilAglIoxx1.1_GCA_905147125.1_num.fa.fai", "rt") as lf:
#with open(alt_fasta_index, "rt") as lf:
    chromLenDict = dict([[s[0],int(s[1])] for s in [l.split() for l in lf]])

long_to_short = chromosome_dict_maker(chromosomes, 2 , '\t', autosomes_only=True)
short_to_long = chromosome_dict_maker(chromosomes, 2 , '\t', autosomes_only=True, short_to_long=True)

print('Creating dataframes of the reference and alternate exon positions with associated gene and mRNA id"s')
#subsetting ref and alt gff to have only gene, mrna and exon positions as well as gene names
ref_exon_gene_pos= ref_annotation.filter_feature_of_type(['seq_id', 'gene','mRNA','exon']).attributes_to_columns()[['seq_id', 'type', 'start', 'end', 'ID', 'Parent', 'strand']]
alt_exon_gene_pos= alt_annotation.filter_feature_of_type(['seq_id','gene','mRNA','exon']).attributes_to_columns()[['seq_id', 'type', 'start', 'end', 'ID', 'Parent', 'strand']]

#getting alt and ref exon only positions dataframe + parent gene and mrna names for each exon 
ref_exon_id=pd.DataFrame(data=exon_gene_mrna_position_finder(ref_exon_gene_pos) ,columns=['chromosome', 'feature_type', 'start', 'end', 'gene_ID', 'mRNA_ID', 'strand'])
alt_exon_id=pd.DataFrame(data=exon_gene_mrna_position_finder(alt_exon_gene_pos) ,columns=['chromosome', 'feature_type', 'start', 'end', 'gene_ID', 'mRNA_ID', 'strand'])

ref_exon_ID_no_duplicates = ref_exon_id.drop_duplicates(subset=['chromosome', 'start', 'end'])
alt_exon_ID_no_duplicates = alt_exon_id.drop_duplicates(subset=['chromosome', 'start', 'end'])

ref_exon_ID_no_duplicates['numeric_chrom'] = ref_exon_ID_no_duplicates['chromosome'].astype(str).apply(long_short_switch)

#getting rid of NA chromosomes that we don't care about
ref_exon_ID_no_duplicates=ref_exon_ID_no_duplicates[ref_exon_ID_no_duplicates['numeric_chrom']!= 'Na']
alt_exon_ID_no_duplicates=alt_exon_ID_no_duplicates[alt_exon_ID_no_duplicates['chromosome']!= 'Na']

#from the gff some of the chromosomes were imported as strings need to be ints
ref_exon_ID_no_duplicates['chromosome'] = ref_exon_ID_no_duplicates['chromosome'].astype(int)


#uncollapsed
total_number_of_unique_exons_ref=len(ref_exon_ID_no_duplicates)
total_number_of_unique_exons_alt=len(alt_exon_ID_no_duplicates)

#collapsed
##needs a version of quality checked exons that can work for the one 
alt_dict = {alt_gene_identification: exons for alt_gene_identification , exons in alt_exon_ID_no_duplicates.groupby('gene_ID')}
alt_dict_collapsed = quaility_check_exons(alt_dict, False)
total_number_of_unique_collapsed_exons_alt = dict_values_to_df(alt_dict_collapsed).sort_values('chromosome').drop('Index', axis = 1)

ref_dict = {ref_gene_identification: exons for ref_gene_identification , exons in ref_exon_ID_no_duplicates.groupby('gene_ID')}
ref_dict_collapsed = quaility_check_exons(ref_dict, False)
total_number_of_unique_collapsed_exons_ref = dict_values_to_df(ref_dict_collapsed).sort_values('chromosome').drop('Index', axis = 1)

#collapsed
total_number_of_unique_exons_collapsed_ref=len(total_number_of_unique_collapsed_exons_ref)
total_number_of_unique_exons_collapsed_alt=len(total_number_of_unique_collapsed_exons_alt)

reference_exons_removed_through_overlapping_etc = total_number_of_unique_exons_ref - total_number_of_unique_exons_collapsed_ref
proportion_reference_exons_removed_through_overlapping_etc = (total_number_of_unique_exons_ref - total_number_of_unique_exons_collapsed_ref)/total_number_of_unique_exons_ref

alternate_exons_removed_through_overlapping_etc = total_number_of_unique_exons_alt - total_number_of_unique_exons_collapsed_alt
proportion_alternate_exons_removed_through_overlapping_etc = (total_number_of_unique_exons_alt - total_number_of_unique_exons_collapsed_alt)/total_number_of_unique_exons_alt

#getting the number of transcripts and exons + names for all genes in the ref and alt assembilies
ref_gene_id_dict = {gene: df for gene, df in ref_exon_ID_no_duplicates.groupby('gene_ID')}    
alt_gene_id_dict = {gene: df for gene, df in alt_exon_ID_no_duplicates.groupby('gene_ID')}

#subetting the above dict to have specific nested dictionary with transcript information etc
#alt_gene_info_dict = gene_transcript_exon_info(alt_gene_id_dict)
#ref_gene_info_dict = gene_transcript_exon_info(ref_gene_id_dict)


print('Determining how many of the reference and alternate exons are alignable')
exon_doest_align_with_alt=0
ref_exon_not_in_aligned_region=0
# Pre-index alt_exon_id by chromosome for faster lookups
alt_exon_id_dict = {chrom : df for chrom, df in alt_exon_ID_no_duplicates.groupby('chromosome')}
ref_exon_id_dict = {chrom: df for chrom, df in ref_exon_ID_no_duplicates.groupby('chromosome')}



#want to add to only take mapping quality of 60
#pysam - when you pull out the read that object has a function called get.aligned.pairs - find the start and end positions - once you have the two numbers there should be an exon in exactly those two position
#gff = 1 based so all exon positions are one based
#bam file = 0 based

all_exons=[]
exons_that_dont_align = []
exons_present_on_aligned_read = []
for chromosome, numeric_chrom in tqdm(short_to_long.items(), total = len(short_to_long)):
    all_reads = list(bamfile.fetch(chromosome))
    for read in all_reads:
        if read_quaility_checks(read, short_to_long, alt_exon_id_dict, 60) == False:
            continue
        ref_exons_present_on_read = ref_exon_id_dict[int(read.reference_name)].loc[(ref_exon_id_dict[int(read.reference_name)]['start'] >= (read.reference_start + 1)) & (ref_exon_id_dict[int(read.reference_name)]['end'] <= (read.reference_end))]
        if len(ref_exons_present_on_read) > 0:
            for exon in ref_exons_present_on_read.itertuples(index=False):
                exons_present_on_aligned_read.append([exon.chromosome, exon.feature_type, exon.start, exon.end, exon.gene_ID, exon.mRNA_ID, exon.strand, exon.numeric_chrom])
        if read.is_forward:
            if read.cigartuples[0][0] != 5:
                if len(read.query_name.split('|')) != 1:
                    alt_exons_present_on_read = alt_exon_id_dict[read.query_name.split('|')[2]].loc[(alt_exon_id_dict[read.query_name.split('|')[2]]['start'] >= (read.query_alignment_start + 1)) & (alt_exon_id_dict[read.query_name.split('|')[2]]['end'] <= (read.query_alignment_end))]
                else:
                    alt_exons_present_on_read = alt_exon_id_dict[read.query_name].loc[(alt_exon_id_dict[read.query_name]['start'] >= (read.query_alignment_start + 1)) & (alt_exon_id_dict[read.query_name]['end'] <= (read.query_alignment_end))]
                pairs_dict={(ref_position+1): (alt_position+1) for alt_position, ref_position in read.get_aligned_pairs(matches_only = False) if alt_position != None if ref_position != None}
            else:
                #if it is soft clipped then we are fine! It is taken care of in the query alignment start bit. but we do need to add hard clipping!
                clipping = read.cigartuples[0][1]
                relative_alt_start = clipping + (read.query_alignment_start+1)
                relative_alt_end = clipping + read.query_alignment_end
                if len(read.query_name.split('|')) != 1:
                    alt_exons_present_on_read = alt_exon_id_dict[read.query_name.split('|')[2]].loc[(alt_exon_id_dict[read.query_name.split('|')[2]]['start'] >= (relative_alt_start)) & (alt_exon_id_dict[read.query_name.split('|')[2]]['end'] <= (relative_alt_end))]
                else:
                    alt_exons_present_on_read = alt_exon_id_dict[read.query_name].loc[(alt_exon_id_dict[read.query_name]['start'] >= (relative_alt_start)) & (alt_exon_id_dict[read.query_name]['end'] <= (relative_alt_end))]
                pairs_dict={(ref_position+1): (clipping + alt_position+1) for alt_position, ref_position in read.get_aligned_pairs(matches_only = False) if alt_position != None if ref_position != None}
        else:
            total_length_query_contig = chromLenDict[read.query_name]
            if read.cigartuples[0][0] != 5:
                #for getting relative positions you need to do 1 based end - 0 based actual positions to end up at the 1 based coordinate rather than 0 based!
                relative_alt_start = total_length_query_contig - (read.query_alignment_end - 1)
                relative_alt_end = total_length_query_contig - read.query_alignment_start
                pairs_dict={(ref_position+1) : ((total_length_query_contig - alt_position)) for alt_position, ref_position in read.get_aligned_pairs(matches_only = False) if alt_position != None and ref_position != None}
            else:
                #if there is hard clipping - you need to remove it using the clipping at the start of the cigar string as it is reversed. So the start is the end and the end is the start relative to the alternate fasta
                if read.cigartuples[1][0] == 5:
                    print('true')
                    clipping=total_length_query_contig  - read.cigartuples[0][1] - read.cigartuples[1][1] 
                    relative_alt_start = clipping - (read.query_alignment_end -1)
                    relative_alt_end = clipping - read.query_alignment_start
                    #make a dictionary of the aligned positions with reference:alt key value format
                    pairs_dict={(ref_position+1) : ((clipping - alt_position)) for alt_position, ref_position in read.get_aligned_pairs(matches_only = False) if alt_position != None and ref_position != None}
                else:
                    clipping=total_length_query_contig  - read.cigartuples[0][1]
                    relative_alt_start = clipping - (read.query_alignment_end -1)
                    relative_alt_end = clipping - read.query_alignment_start
                    #make a dictionary of the aligned positions with reference:alt key value format
                    pairs_dict={(ref_position+1) : ((clipping - alt_position)) for alt_position, ref_position in read.get_aligned_pairs(matches_only = False) if alt_position != None and ref_position != None}
            #pairs_dict={(ref_position+1) : ((total_length_query_contig - alt_position)) for alt_position, ref_position in read.get_aligned_pairs(matches_only = False) if alt_position != None and ref_position != None}
            #alt_exons_present_on_read = alt_exon_id_dict[read.query_name].loc[((alt_exon_id_dict[read.query_name]['start']) >= (relative_alt_start))& (alt_exon_id_dict[read.query_name]['end'] <= (relative_alt_end))]
            if len(read.query_name.split('|')) != 1:
                alt_exons_present_on_read = alt_exon_id_dict[read.query_name.split('|')[2]].loc[((alt_exon_id_dict[read.query_name.split('|')[2]]['start']) >= (relative_alt_start))& (alt_exon_id_dict[read.query_name.split('|')[2]]['end'] <= (relative_alt_end))]
            else:
                alt_exons_present_on_read = alt_exon_id_dict[read.query_name].loc[((alt_exon_id_dict[read.query_name]['start']) >= (relative_alt_start))& (alt_exon_id_dict[read.query_name]['end'] <= (relative_alt_end))]
        alt_exons_on_read_group_by_gene = {geneID : info for geneID, info in alt_exons_present_on_read.groupby('gene_ID')}
        ref_exons_present_on_read_group_by_gene = {geneID : info for geneID, info in ref_exons_present_on_read.groupby('gene_ID')}
        for key, exons_df in ref_exons_present_on_read_group_by_gene.items():
            #is the alternate version of the gene part of the exons present on the aligned contig?
            if key in alt_exons_on_read_group_by_gene.keys():
                alternate_version_of_gene = alt_exons_on_read_group_by_gene[key]
                #for each exon determining if it has a reciprocal position in the alternate
                for values in exons_df.itertuples(index = False):
                    if read.is_forward:
                        #are the reference exons on this read?
                        if values.start in pairs_dict and values.end in pairs_dict:
                            #if it is what are the reciprocal positions in the alternate?
                            alt_start, alt_end = pairs_dict[values.start] , pairs_dict[values.end]
                        else:
                            ref_exon_not_in_aligned_region += 1
                            continue
                    else:
                        if values.start in pairs_dict and values.end in pairs_dict:
                            alt_end , alt_start = pairs_dict[values.start] , pairs_dict[values.end]
                        else:
                            ref_exon_not_in_aligned_region += 1
                            continue
                    #is there an alternate exon aligned perfectly?
                    final_one = alternate_version_of_gene.loc[((alternate_version_of_gene.start) == alt_start) & ((alternate_version_of_gene.end) == alt_end)]
                    if len(final_one) > 0:
                        for values2 in final_one.itertuples(index = False):
                            all_exons.append([values.chromosome, values.start, values.end, values.gene_ID, values.mRNA_ID, values.strand, values.numeric_chrom, values2.chromosome, values2.feature_type, values2.start, values2.end, values2.gene_ID, values2.mRNA_ID, values2.strand, 'both', 'NA', read.is_forward])
                    elif len(final_one) == 0:
                        #determine if the start aligns but not the end
                        final_one = alternate_version_of_gene.loc[((alternate_version_of_gene.start) == alt_start)]
                        for values2 in final_one.itertuples(index = False):
                            all_exons.append([values.chromosome, values.start, values.end, values.gene_ID, values.mRNA_ID, values.strand, values.numeric_chrom, values2.chromosome, values2.feature_type, values2.start, values2.end, values2.gene_ID, values2.mRNA_ID, values2.strand, 'start', alt_end, read.is_forward])
                    elif len(final_one) == 0:
                        #determine if the end aligns but not the start
                        final_one = alternate_version_of_gene.loc[((alternate_version_of_gene.end) == alt_end)]
                        for values2 in final_one.itertuples(index = False):
                            all_exons.append([values.chromosome, values.start, values.end, values.gene_ID, values.mRNA_ID, values.strand, values.numeric_chrom, values2.chromosome, values2.feature_type, values2.start, values2.end, values2.gene_ID, values2.mRNA_ID, values2.strand, 'end', alt_start, read.is_forward])
                    if len(final_one) == 0:
                        #determine if neither of the ends align
                        exons_that_dont_align.append([values.chromosome, values.start, values.end, values.gene_ID, values.mRNA_ID, values.strand, values.numeric_chrom, read.is_forward])
                        #find the gene in the alternate
                        alternate_gene_version = alt_gene_id_dict[values.gene_ID]
                        #find the reciprocal positions
                        #pairs_dict[values.start] , pairs_dict[values.end]

                    
#start=87583
#end=88040
#gene_ID = 'gene:ENSBPTG00005003567'

#exons that aren't present in aligned regions of the reference and the alterante
exons_that_are_in_aligned_regions_df = pd.DataFrame(exons_present_on_aligned_read, columns=['ref_chrom', 'feature_type', 'ref_start', 'ref_end', 'ref_gene_id', 'ref_transcript_id', 'ref_strand', 'ref_chrom_num'])

#exons that don't align
exons_that_dont_align_df = pd.DataFrame(exons_that_dont_align, columns=['ref_chrom', 'ref_start', 'ref_end', 'ref_gene_id', 'ref_transcript_id', 'ref_strand', 'ref_chrom_num', 'contig_alignment_orientation'])

#alt_gene_id_dict['gene:ENSBPTG00005003567']
#ref_gene_id_dict['gene:ENSBPTG00005003567'] 

#pull out any snp differences in that region -> could just throw away at that point
#might not be one for in aligining 

all_exons_df = pd.DataFrame(all_exons, columns=['ref_chrom', 'ref_start', 'ref_end', 'ref_gene_id', 'ref_transcript_id', 'ref_strand', 'ref_chrom_num', 'alt_chrom', 'alt_feature', 'alt_start', 'alt_end', 'alt_gene_id', 'alt_transcript_id', 'alt_strand', 'aligned', 'non_aligned_end', 'contig_alignment_orientation'])

not_both_aligned_exons_df=all_exons_df[all_exons_df['aligned'] != 'both']

not_both_aligned_exons_df = not_both_aligned_exons_df[['ref_chrom', 'ref_start', 'ref_end', 'alt_start', 'alt_end', 'aligned', 'non_aligned_end']]

all_exons_df=all_exons_df[all_exons_df['aligned'] == 'both']
#remvoing the mutliple versions from transcripts
#i think there is an argument for removing all versions of a exon that is being aligned to multiple times as we could never determine which gene it is coming from? 
all_exons_df = all_exons_df.drop_duplicates(subset=['ref_gene_id','ref_start', 'ref_end'], keep=False)
all_exons_df = all_exons_df.drop_duplicates(subset=['alt_gene_id','alt_start', 'alt_end'], keep=False)


print((len(all_exons_df)/total_number_of_unique_exons_ref)*100)
print((len(all_exons_df)/total_number_of_unique_exons_alt)*100)

#added in the filter out any nested exons from the REF prespective! 
#this is all the ref perspective 

#For ref genes filter the exons so we have only one alternate contig associated with it and also for nested exons we only end up with the longest one!:

#Removing nested exons for the alternate after already doing it for the reference - makes sure there are no nested - does it also only keep the longest ones?
## only concern is that what if the order of doing it means i miss something 
aligned_exons_alt_focus_dict = {alt_gene_identification: exons for alt_gene_identification , exons in all_exons_df.groupby('alt_gene_id')}
alternate_gene_exons = quaility_check_exons(aligned_exons_alt_focus_dict, 'alt')
alt_df = dict_values_to_df(alternate_gene_exons).sort_values('alt_chrom').drop('Index', axis = 1)

aligned_exons_ref_focus_dict = {ref_gene_identification: exons for ref_gene_identification , exons in alt_df.groupby('ref_gene_id')}
reference_gene_exons = quaility_check_exons(aligned_exons_ref_focus_dict, 'ref')
ref_df = dict_values_to_df(reference_gene_exons).sort_values('ref_chrom').drop('Index', axis = 1)


#so at this point we have removed nested indels that are nested in both ref and alt and we have only kept the longest version of each indel - it has to be the longest in both the ref and alt or it prints a special message.
#so we essentially for each ref gene we have all the genomic regions that are ex-onic that are aligned to exonic tracts in the alt
##distance_from_actual_exon_start_or_end=[]
##proposed_length_of_trying_to_align_exon=[]

##for rows in not_both_aligned_exons_df.itertuples(index=False):
##    if rows.aligned == 'end':
##        distance_from_actual_exon_start_or_end.append(rows.alt_start - rows.non_aligned_end)
##        proposed_length_of_trying_to_align_exon.append(rows.alt_end  - rows.non_aligned_end)
##    elif rows.aligned == 'start':
##        distance_from_actual_exon_start_or_end.append(rows.alt_end - rows.non_aligned_end)
##        proposed_length_of_trying_to_align_exon.append(rows.non_aligned_end - rows.alt_start)
       
##not_both_aligned_exons_df['distance_from_actual_exon_start_or_end'] = distance_from_actual_exon_start_or_end
##not_both_aligned_exons_df['proposed_length_of_trying_to_align_exon'] = proposed_length_of_trying_to_align_exon

#checking that there is no bias for exons to be bias to being shorter or longer in the alternate 
#not_both_aligned_exons_df_subsetted = not_both_aligned_exons_df[not_both_aligned_exons_df['distance_from_actual_exon_start_or_end'] >= -100]
#not_both_aligned_exons_df_subsetted = not_both_aligned_exons_df_subsetted[not_both_aligned_exons_df_subsetted['distance_from_actual_exon_start_or_end'] <= 100]


#stats
average_length_of_unalignable_exon=np.mean(exons_that_dont_align_df['ref_end'] - exons_that_dont_align_df['ref_start'])
average_length_of_alignable_exon = np.mean(ref_df['ref_end'] - ref_df['ref_start'])
total_number_of_aligned_exon_bp=sum(ref_df['ref_end'] - ref_df['ref_start'])
what_proportion_of_possible_exonic_bp_are_aligned=total_number_of_aligned_exon_bp/(sum(total_number_of_unique_collapsed_exons_ref['end'] - total_number_of_unique_collapsed_exons_ref['start']))

proportion_of_exons_in_unaligned_regions_of_the_genome = (total_number_of_unique_exons_collapsed_ref - len(exons_that_are_in_aligned_regions_df))/total_number_of_unique_exons_collapsed_ref

proportion_of_exons_that_didnt_align_on_both_ends = len(not_both_aligned_exons_df)/total_number_of_unique_exons_collapsed_ref
proportion_of_exons_that_didnt_align_at_all = len(exons_that_dont_align_df)/total_number_of_unique_exons_collapsed_ref
total_proportion_of_reference_exons_that_are_aligned = (len(ref_df)/total_number_of_unique_exons_collapsed_ref)
total_proportion_of_alternate_exons_that_are_aligned = (len(ref_df)/total_number_of_unique_exons_collapsed_alt)
total_proportion_of_nested_exons_collated = (len(all_exons_df)- len(ref_df))/total_number_of_unique_exons_collapsed_ref

Total_assembly_info=[total_proportion_of_reference_exons_that_are_aligned, total_proportion_of_alternate_exons_that_are_aligned,total_proportion_of_nested_exons_collated, what_proportion_of_possible_exonic_bp_are_aligned, proportion_of_exons_in_unaligned_regions_of_the_genome, proportion_of_exons_that_didnt_align_at_all, proportion_of_exons_that_didnt_align_on_both_ends]
Total_assembly_info_df = pd.DataFrame(Total_assembly_info)
Total_assembly_info_df_t = Total_assembly_info_df.transpose()
Total_assembly_info_df_t.columns = ['Proportion_of_reference_exons_aligned', 'Proportion_of_alternate_exons_aligned', 'Proportion_of_nested_exons_collapsed', 'What_proportion_of_reference_exonic_bp_are_aligned','Proportion_of_exons_in_unalignable_regions_of_the_genome' ,'Proportion_of_exons_that_didnt_align_between_haplotypes', 'Proportion_of_exons_that_aligned_one_end']


print(f'The total percentage of reference exons that are alignable is {total_proportion_of_reference_exons_that_are_aligned}%')
print(f'The total percentage of alternate exons that are alignable is {total_proportion_of_alternate_exons_that_are_aligned}%')

#ref_df.to_csv('/home/s1929681/Desktop/filtered_aligned_exons.csv' , index = False)
Total_assembly_info_df_t.to_csv(args.exon_alignment_stats, index = False)

ref_df.to_csv(args.aligned_exons_output, index = False)

#5% on either side
#ploy that % of the exons is the actual missing chunks 
#check the distribution of evens and positives 

## Plotting side plot
#import seaborn as sns
#import matplotlib.pyplot as plt
#unalinged_ends = all_exons_df[all_exons_df['aligned'] == 'start']
#unalinged_ends['difference_between_ends'] = unalinged_ends['non-alinged_end'] - unalinged_ends['alt_end']
#unaligned_start = all_exons_df[all_exons_df['aligned'] == 'end']
#unaligned_start['difference_between_starts'] = unaligned_start['non-alinged_end'] - unaligned_start['alt_start']

#all_exons_df_subsetted = all_exons_df[all_exons_df['alt_chrom'] == 'CAKNFC010000229.1']
#all_exons_df_subsetted_all_aligned = all_exons_df_subsetted[all_exons_df_subsetted['aligned'] == 'both']
#sns.relplot(data=all_exons_df_subsetted, x ='ref_start', y='ref_end',hue='ref_gene_id' )

#sns.relplot(data=all_exons_df_subsetted, x ='alt_start', y='alt_end',hue='ref_gene_id' )

#sns.relplot(data=all_exons_df_subsetted, x ='ref_start', y='alt_start',hue='ref_gene_id' )

#sns.relplot(data=all_exons_df_subsetted, x ='ref_start', y='alt_start',hue='alt_chrom' )

#sns.relplot(data=all_exons_df_subsetted, x ='ref_gene_id', y='alt_start',hue='ref_gene_id' )

#sns.relplot(data=all_exons_df_subsetted, x ='ref_gene_id', y='ref_start',hue='ref_gene_id' )
#plt.show()
