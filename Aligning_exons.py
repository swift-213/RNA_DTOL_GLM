#okay so the aim is to take the exon positions in both the ref and the alt and determine if they are actually aligned in the haplotypes and if they are aligned if the region is classes as a exon in the alt genome gff!
#So we will need the ref and alt bamfiles to ensure they are aligned as well as the ref and alt gffs!
import pysam
import gffpandas.gffpandas as gffpd
import pandas as pd
from tqdm import tqdm
from genfun import line_splitter
import numpy as np

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
        meta_info = [no_transcripts, number_exons, mRNA_ids,exon_positions, chromosome]
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
    elif reads.query_name.split('|')[2] not in alt_chromosome_name_dict:
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

def quaility_check_exons(gene_exon_dict, reference_or_alternate):
    No_duplicate_start_end_positions={}
    all_gene_ids_with_the_aligned_ids = {}
    thrown_out_ids = set()
    string_gene_id = str(reference_or_alternate) + '_gene_id'
    chromosome = str(reference_or_alternate) + '_chrom'
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
            if not any((exons.ref_start < row.ref_start) & (exons.ref_end > row.ref_end) & (exons.index != index)):
                non_nested_exons.append(index)
            if not any((exons.alt_start < row.alt_start) & (exons.alt_end > row.alt_end) & (exons.index != index)):
                non_nested_exons.append(index)
        #exon_counter_dict = counter(non_nested_exons)
        #all_non_nested_exons = removing_non_duplicated_values(non_nested_exons, exon_counter_dict)
        exons = (exons.loc[list(set(non_nested_exons))]).sort_index()
        #here we are trying to remove exons that have the same start or end position and taking only the longest one
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
    return No_duplicate_start_end_positions

def dict_values_to_df(dictionary):
    values_as_list=[]
    for key, value in dictionary.items():
        for row in value.itertuples():
            values_as_list.append(row)
    return pd.DataFrame(values_as_list)


print('Loading in files')
chromosomes = ("/home/s1929681/One_drive_file_copies_25_01_15/RA_work/Indel_Project/data_set_chroms/iyBomPrat1.1_chromosomes.txt")
#Reading in ref bam
bamfile = pysam.AlignmentFile("/media/s1929681/Seagate/Frankie_DTOL_lep_project/outputs/samtools/iyBomPrat1.1_alignment.sort.bam", "rb")
#Reading in ref gff
ref_annotation = gffpd.read_gff3('/media/s1929681/Seagate/Frankie_DTOL_lep_project/phase_2_gff/ref_gff/iyBomPrat1.1_GCA_930367275.1.gff3')
#Reading in alt gff
alt_annotation = gffpd.read_gff3('/media/s1929681/Seagate/Frankie_DTOL_lep_project/phase_2_gff/alt_gff/iyBomPrat1.1_alt.gff3')

#Getting the start and end positions for each chromosome in order to do the array approach for SNP and monomorphic sites 
with open("/media/s1929681/Seagate/Frankie_DTOL_lep_project/phase_two_fasta/iyBomPrat1.1_GCA_930367225.1.fa.fai", "rt") as lf:
#with open(fasta_fai, "rt") as lf:
    chromLenDict = dict([[s[0],int(s[1])] for s in [l.split() for l in lf]])

chrom_key=[]
line_splitter(chromosomes, 2, chrom_key, '\t')
chrom_df=pd.DataFrame(chrom_key, columns=['Long', 'short'])

long_to_short = chrom_df.set_index("Long")["short"].to_dict()
short_to_long = chrom_df.set_index("short")["Long"].to_dict()

chrom_df['Long']=chrom_df['Long'].astype("string")
chrom_df['short']=chrom_df['short'].astype("string")

print('Creating dataframes of the reference and alternate exon positions with associated gene and mRNA id"s')
#subsetting ref and alt gff to have only gene, mrna and exon positions as well as gene names
ref_exon_gene_pos= ref_annotation.filter_feature_of_type(['seq_id', 'gene','mRNA','exon']).attributes_to_columns()[['seq_id', 'type', 'start', 'end', 'ID', 'Parent', 'strand']]
alt_exon_gene_pos= alt_annotation.filter_feature_of_type(['seq_id','gene','mRNA','exon']).attributes_to_columns()[['seq_id', 'type', 'start', 'end', 'ID', 'Parent', 'strand']]

#getting alt and ref exon only positions dataframe + parent gene and mrna names for each exon 
ref_exon_id=pd.DataFrame(data=exon_gene_mrna_position_finder(ref_exon_gene_pos) ,columns=['chromosome', 'feature_type', 'start', 'end', 'gene_ID', 'mRNA_ID', 'strand'])
alt_exon_id=pd.DataFrame(data=exon_gene_mrna_position_finder(alt_exon_gene_pos) ,columns=['chromosome', 'feature_type', 'start', 'end', 'gene_ID', 'mRNA_ID', 'strand'])
            
ref_exon_id['numeric_chrom'] = ref_exon_id['chromosome'].astype(str).apply(long_short_switch)

#getting the number of transcripts and exons + names for all genes in the ref and alt assembilies
ref_gene_id_dict = {gene: df for gene, df in ref_exon_id.groupby('gene_ID')}    
alt_gene_id_dict = {gene: df for gene, df in alt_exon_id.groupby('gene_ID')}

alt_gene_info_dict = gene_transcript_exon_info(alt_gene_id_dict)
ref_gene_info_dict = gene_transcript_exon_info(ref_gene_id_dict)

#output the ones that fail as well 
#pull out any snp differences in that region -> could just throw away at that point
#might not be one for in aligining 

print('Determining how many of the reference and alternate exons are alignable')
exon_doest_align_with_alt=0
ref_exon_not_in_aligned_region=0
more_than_one_exon_aligning_to_one_read=0
aligned_exons=[]
mutliple_exon_alignment_hits=[]
#all_ref_exons_looked_through = set()
#all_alt_exons_looked_through = set()
# Pre-index alt_exon_id by chromosome for faster lookups
alt_exon_id_dict = {chrom : df for chrom, df in alt_exon_id.groupby('chromosome')}
ref_exon_id_dict = {chrom: df for chrom, df in ref_exon_id.groupby('numeric_chrom')}

#want to add to only take mapping quality of 60
#pysam - when you pull out the read that object has a function called get.aligned.pairs - find the start and end positions - once you have the two numbers there should be an exon in exactly those two positions
all_reads = list(bamfile.fetch())
for read in tqdm(all_reads, total = len(all_reads)):
    if read_quaility_checks(read, long_to_short, alt_exon_id_dict, 60) == False:
        continue
    ref_exons_present_on_read = ref_exon_id_dict[read.reference_name].loc[(ref_exon_id_dict[read.reference_name]['start'] >= (read.reference_start + 1)) & (ref_exon_id_dict[read.reference_name]['end'] <= (read.reference_end))]
    if read.is_forward:
        alt_exons_present_on_read = alt_exon_id_dict[read.query_name.split('|')[2]].loc[(alt_exon_id_dict[read.query_name.split('|')[2]]['start'] >= (read.query_alignment_start + 1)) & (alt_exon_id_dict[read.query_name.split('|')[2]]['end'] <= (read.query_alignment_end))]
        pairs_dict={(ref_position+1): (alt_position+1) for alt_position, ref_position in read.get_aligned_pairs(matches_only = False) if alt_position != None if ref_position != None}
    else:
        total_length_query_contig = chromLenDict[read.query_name]
        if read.cigartuples[-1][0] != 4 and read.cigartuples[-1][0] != 5:
            #for getting relative positions you need to do 1 based end - 0 based actual positions to end up at the 1 based coordinate rather than 0 based!
            relative_alt_start = total_length_query_contig - (read.query_alignment_end - 1)
            relative_alt_end = total_length_query_contig - read.query_alignment_start
        else:
            #what if there is hard and soft clipping? 
            relative_alt_start = (total_length_query_contig - read.cigartuples[-1][1]) - (read.query_alignment_end -1)
            relative_alt_end = (total_length_query_contig - read.cigartuples[-1][1]) - read.query_alignment_start
        alt_exons_present_on_read = alt_exon_id_dict[read.query_name.split('|')[2]].loc[((alt_exon_id_dict[read.query_name.split('|')[2]]['start']) >= (relative_alt_start))& (alt_exon_id_dict[read.query_name.split('|')[2]]['end'] <= (relative_alt_end))]
        #make a dictionary of the aligned positions with reference:alt key value format
        pairs_dict={(ref_position+1) : ((total_length_query_contig - alt_position)) for alt_position, ref_position in read.get_aligned_pairs(matches_only = False) if alt_position != None and ref_position != None}
    for values in ref_exons_present_on_read.itertuples(index = False):
        if read.is_forward:
            if values.start in pairs_dict and values.end in pairs_dict:
                alt_start, alt_end = pairs_dict[values.start] , pairs_dict[values.end]
            else:
                ref_exon_not_in_aligned_region += 1
        else:
            if values.start in pairs_dict and values.end in pairs_dict:
                alt_end , alt_start = pairs_dict[values.start] , pairs_dict[values.end]
            else:
                ref_exon_not_in_aligned_region += 1
        final_one = alt_exons_present_on_read.loc[((alt_exons_present_on_read.start) == alt_start) & ((alt_exons_present_on_read.end) == alt_end)]
        if len(final_one) > 0:
            if len(final_one) > 1:
                more_than_one_exon_aligning_to_one_read += 1
                appending_the_exon_df(final_one, mutliple_exon_alignment_hits)
            else:
                appending_the_exon_df(final_one, aligned_exons) 
        else:
            exon_doest_align_with_alt +=1

#turning the lists into dataframes
exon_double_bed_df = pd.DataFrame(aligned_exons, columns=['ref_chrom', 'ref_chrom_num', 'ref_feature', 'ref_start', 'ref_end', 'ref_gene_id', 'ref_transcript_id', 'ref_strand', 'alt_chrom', 'alt_feature', 'alt_start', 'alt_end', 'alt_gene_id', 'alt_transcript_id', 'alt_read_forward', 'alt_strand'])
multiple_hits_bby = pd.DataFrame(mutliple_exon_alignment_hits, columns=['ref_chrom', 'ref_chrom_num', 'ref_feature', 'ref_start', 'ref_end', 'ref_gene_id', 'ref_transcript_id', 'ref_strand', 'alt_chrom', 'alt_feature', 'alt_start', 'alt_end', 'alt_gene_id', 'alt_transcript_id', 'alt_read_forward', 'alt_strand'])

#join all the df's together! + removing all the duplicated rows from different transcripts 
all_single_copy_exon_df = pd.concat([exon_double_bed_df, multiple_hits_bby], ignore_index = True)
duplicates_ignoring_alt_and_ref_transcripts = all_single_copy_exon_df.duplicated(subset = all_single_copy_exon_df.columns.difference(['alt_transcript_id', 'ref_transcript_id']))
all_single_copy_exon_df = all_single_copy_exon_df[~duplicates_ignoring_alt_and_ref_transcripts]

all_single_copy_exon_df = all_single_copy_exon_df.sort_values(['ref_chrom', 'ref_start'], ascending=[True, True], ignore_index=True)

total_number_of_reference_exons_that_are_aligned = (len(all_single_copy_exon_df)/len(ref_exon_gene_pos))*100
total_number_of_alternate_exons_that_are_aligned = (len(all_single_copy_exon_df)/len(alt_exon_gene_pos))*100
print(f'The total percentage of reference exons that are alignable is {total_number_of_reference_exons_that_are_aligned}%')
print(f'The total percentage of alternate exons that are alignable is {total_number_of_alternate_exons_that_are_aligned}%')
#added in the filter out any nested exons from the REF prespective! 
#this is all the ref perspective 

#def counter(list):
#    count_dict={}
#    for i in list:
#        count_dict[i] = list.count(i)
#    return count_dict

#def removing_non_duplicated_values(list_for_checking, counter_dict):
#    tested = set()
#    for value in list_for_checking:
#        if value not in tested:
            #if counter_dict[value] != 1:
#            tested.add(value)
#    return list(tested)



#so at this point we have removed nested indels that are nested in both ref and alt and we have only kept the longest version of each indel - it has to be the longest in both the ref and alt or it prints a special message.
#so we essentially for each ref gene we have all the genomic regions that are ex-onic that are aligned to exonic tracts in the alt
#For ref genes filter the exons so we have only one alternate contig associated with it and also for nested exons we only end up with the longest one!:
aligned_exons_ref_focus_dict = {ref_gene_identification: exons for ref_gene_identification , exons in all_single_copy_exon_df.groupby('ref_gene_id')}
aligned_exons_alt_focus_dict = {alt_gene_identification: exons for alt_gene_identification , exons in all_single_copy_exon_df.groupby('alt_gene_id')}


reference_gene_exons = quaility_check_exons(aligned_exons_ref_focus_dict, 'alt')
len(reference_gene_exons)


alternate_gene_exons = quaility_check_exons(aligned_exons_alt_focus_dict, 'ref')
len(alternate_gene_exons)



alt_df = dict_values_to_df(alternate_gene_exons).sort_values('alt_chrom')
ref_df = dict_values_to_df(reference_gene_exons).sort_values('ref_chrom').drop('Index', axis = 1)

#I need to work out how many exons we are loosing at each stage
    #going to compare to the orginal number of exons in the gff files! Will say for alt exons we retained X amount and for ref exons we retained X amount but will do this once we have the final count!

orginal_number_of_exons_in_ref = len(ref_exon_gene_pos)
final_number_of_exons_that_are_aligned_and_filtered = len(ref_df)

total_percentage_of_reference_exons_kept = (final_number_of_exons_that_are_aligned_and_filtered/orginal_number_of_exons_in_ref)*100
total_percentage_of_alternate_exons_kept = (final_number_of_exons_that_are_aligned_and_filtered/len(alt_exon_gene_pos))*100
print(f'After filtering the total percentage of alignable reference exons is {total_percentage_of_reference_exons_kept}%. This is {final_number_of_exons_that_are_aligned_and_filtered} exons in total')
print(f'After filtering the total percentage of alignable reference exons is {total_percentage_of_alternate_exons_kept}%. This is {final_number_of_exons_that_are_aligned_and_filtered} exons in total')

ref_df.to_csv('/home/s1929681/Desktop/filtered_aligned_exons.csv' , index = False)
