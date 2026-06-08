import gffpandas.gffpandas as gffpd
import pandas as pd
from argparse import ArgumentParser


def exon_gene_mrna_position_finder(input_df):
    exon_id=[]
    gene_ID = None  
    mRNA_ID = None
    for row, values in input_df.iterrows():
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


def sort_by_gene_output(input_gene_id_dict):
    per_gene_output=[]
    for key, value in input_gene_id_dict.items():
        number_of_transcripts=value['mRNA_ID'].nunique()
        number_of_exons=len(value)
        number_of_unique_exon_starts=value['start'].nunique()
        number_of_unique_exon_ends=value['end'].nunique()
        number_of_unique_exons=len(value.drop_duplicates(subset=['chromosome', 'start', 'end']))
        strand=value['strand'].nunique()
        per_gene_output.append([value['chromosome'].iloc[0], key, min(value['start']), max(value['end']), number_of_unique_exons, number_of_transcripts, number_of_unique_exon_starts, number_of_unique_exon_ends,number_of_exons, strand])
    df_per_gene_output = pd.DataFrame(per_gene_output, columns=['Chromosome', 'Gene_ID', 'First_exon_pos', 'Last_exon_end', 'Number_of_unique_exons', 'Number_of_transcripts', 'Number_of_unique_exon_starts', 'Number_of_unique_exon_ends', 'Total_exons_associated', 'Strand'])
    return df_per_gene_output

parser = ArgumentParser()
#inputs
parser.add_argument("--reference_gff", help="reference gff file unzipped", required=True)
parser.add_argument("--alternate_gff", help="alternate gff file unzipped", required=True)


#outputs
parser.add_argument("--reference_gff_information", action='store', help="csv file output of every reference gene with gff metadata", required=True)
parser.add_argument("--alternate_gff_information", action='store', help="csv file output of every alternate gene with gff metadata", required=True)

args = parser.parse_args()

##reading in files
ref_annotation = gffpd.read_gff3(args.reference_gff)
alt_annotation = gffpd.read_gff3(args.alternate_gff)

##ref_annotation = gffpd.read_gff3('/mnt/loki/martin/frankie/RNA_seq_glm/raw_files/gff_files/icAbaPara2.1_GCA_964197645.1.gff3')
#Reading in alt gff
##alt_annotation = gffpd.read_gff3('/mnt/loki/martin/frankie/RNA_seq_glm/raw_files/gff_files/icAbaPara2.1_GCA_964197635.1_liftover.gff3')


print('Creating dataframes of the reference and alternate exon positions with associated gene and mRNA id"s')
#subsetting ref and alt gff to have only gene, mrna and exon positions as well as gene names
ref_exon_gene_pos= ref_annotation.filter_feature_of_type(['seq_id', 'gene','mRNA','exon']).attributes_to_columns()[['seq_id', 'type', 'start', 'end', 'ID', 'Parent', 'strand']]
alt_exon_gene_pos= alt_annotation.filter_feature_of_type(['seq_id','gene','mRNA','exon']).attributes_to_columns()[['seq_id', 'type', 'start', 'end', 'ID', 'Parent', 'strand']]

#getting alt and ref exon only positions dataframe + parent gene and mrna names for each exon 
ref_exon_id=pd.DataFrame(data=exon_gene_mrna_position_finder(ref_exon_gene_pos) ,columns=['chromosome', 'feature_type', 'start', 'end', 'gene_ID', 'mRNA_ID', 'strand'])
alt_exon_id=pd.DataFrame(data=exon_gene_mrna_position_finder(alt_exon_gene_pos) ,columns=['chromosome', 'feature_type', 'start', 'end', 'gene_ID', 'mRNA_ID', 'strand'])


ref_gene_dict = {ref_gene_identification: exons for ref_gene_identification , exons in ref_exon_id.groupby('gene_ID')}
alt_gene_dict = {alt_gene_identification: exons for alt_gene_identification , exons in alt_exon_id.groupby('gene_ID')}

ref_gene_information = sort_by_gene_output(ref_gene_dict)
alt_gene_information = sort_by_gene_output(alt_gene_dict)

ref_gene_information.to_csv(args.reference_gff_information, index=False)
alt_gene_information.to_csv(args.alternate_gff_information, index=False)

