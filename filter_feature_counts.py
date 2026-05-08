import pandas as pd
from general.genfun import chromosome_dict_maker
from argparse import ArgumentParser

parser = ArgumentParser()
#inputs
parser.add_argument("-chromosomes", "--chromosome_name_files", help="File containing the chromosome names for the reference genome", required=True)
parser.add_argument("-feature_counts", "--feature_counts", help="The orginal unfiltered feature counts output", required=True)

#outputs
parser.add_argument("-filtered_feature_counts", "--filtered_feature_counts", action='store', help="The filtered for autosomes only feature counts", required=True)

args = parser.parse_args()

##chromosomes = ("/home/s1929681/git_directories/ac3/Chromosome_name_conversion_files/ilGriApri1.1_chromosomes.txt")
##feature_counts = pd.read_csv('/media/s1929681/Seagate_B/Frankie_DTOL_lep_project/RNA_DTOL_GLM_OUTPUTS/5_feature_counts/ilGriApri1.1_GCA_916610205.1_pcrdups_aligned_exon_annotated_featurecounts.txt', delimiter='\t', skiprows=[0])

chromosomes = (args.chromosome_name_files)
feature_counts = pd.read_csv(args.feature_counts, delimiter='\t', skiprows=[0])


chromosomes_to_keep = chromosome_dict_maker(chromosomes, 2, '\t', autosomes_only=True, short_to_long=True)

filtered_df = feature_counts[
    feature_counts["Chr"].apply(
        lambda x: any(chromosome in chromosomes_to_keep.keys() for chromosome in x.split(";"))
    )
]

filtered_df.to_csv(args.filtered_feature_counts, index = False, sep='\t')