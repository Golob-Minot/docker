#!/usr/bin/env python

# Goal: Create a more human-readable matrix of results
# Here: Each specimen gets a COLUMN (optionally with human-curated labels for the specimens below)
#       Each classified TAXON gets a ROW (Please note: sequence variants colliding into the same taxon are ADDED. This is NOT IDEAL, but reasonable for a quick look.)


import pandas as pd
import csv
import argparse
import sys

def get_best_name_rank(tax_row):
    rank_order = list(tax_row.index)
    name = None
    while rank_order and name is None:
        rank = rank_order.pop()
        if not pd.isnull(tax_row[rank]):
            name = tax_row[rank]
    if rank == 'Species':
        name = tax_row['Genus']+" "+name
    return name, rank        

def main():
    args_parser = argparse.ArgumentParser()
    
    args_parser.add_argument('--seqtable','-s', help="Sequence table from dada2, in CSV format", required =True, type=file)
    args_parser.add_argument('--taxonomy','-t', help="taxonomy for sequence variants in CSV format as from dada2", required =True, type=file)
    args_parser.add_argument('--tallies_wide', '-T', help="Where to put the tallies-wide format output", type=argparse.FileType('w'), required=True)
    
    # optional
    args_parser.add_argument('--fractions_wide', '-F', help="(OPTIONAL) Where to put the fractions-wide format output", type=argparse.FileType('w'))
    args_parser.add_argument('--labels','-L', help="(optional) Specimens -> Human-annotated labels", type=file)
    
    args = args_parser.parse_args()
    
    # Load the sequence table
    seqtab = pd.read_csv(args.seqtable, index_col=0)
    # Load the taxonomy
    taxonomy = pd.read_csv(args.taxonomy, index_col=0)
    # extract the rank order
    rank_order = taxonomy.columns.tolist()
    
    # seqtab columns are sequence variants. taxonomy row names are the same sequence variants
    
    # first step is to get the best identification for each sequence variant
    best_name_rank = pd.DataFrame()
    best_name_rank['name_rank'] = taxonomy.apply(get_best_name_rank,axis=1)
    best_name_rank['tax_name'] = best_name_rank.name_rank.apply(lambda nr: nr[0])
    best_name_rank['tax_rank'] = best_name_rank.name_rank.apply(lambda nr: nr[1])
    
    # Taxtable is a bit crude. It lumps together everything classified to the same name. Let's do that using dictionary comprehension, where the key is a unique tax_name and value is a nested dict of rank and sequence variants 
    
    tax_names = {tax_name: {
        'tax_rank': t_block.tax_rank.unique()[0],
        'sv':       t_block.index,
    } for tax_name, t_block in best_name_rank.groupby('tax_name')}
    
    
    # Create the tallies-wide dataframe
    tw_df = pd.DataFrame()
    
    for tax_name in tax_names:
        tw_df[tax_name] = seqtab[tax_names[tax_name]['sv']].apply(lambda r: r.sum(), axis=1)
        
    # Add the ranks
    tw_df.loc['rank'] = tw_df.apply(lambda c: tax_names[c.name]['tax_rank'])
    
    # Add the labels (if we are given labels):
    if args.labels:
        labels_df = pd.read_csv(args.labels)
        tw_df = pd.merge(tw_df,labels_df, how='left', left_index=True, right_on='specimen')
        # reset the index
        tw_df.index = tw_df.specimen
        # Get rid of the duplicated specimen column added by the join
        tw_df.drop('specimen', axis=1, inplace=True)
        
        # clean up missing labels
        tw_df['label'].fillna("", inplace=True)
        
        # Sort the columns, putting rank to the far left, then sorting the by rank and name
        taxon_cols = [c for c in tw_df.columns.tolist() if c != 'label']
        # Lookup the rank order for each taxon column, use zip to make a tuple of (rank_order_i, taxon_column_name)
        rank_taxon = zip([rank_order.index(tw_df[c]['rank']) for c in taxon_cols],taxon_cols)
        # Sort the tuples, using sorted that is smart enough to use the tuple properly (first sort by the rank order, then by name)
        # return the column name of teh rest
        sorted_taxon_cols = [rt[1] for rt in sorted(rank_taxon)]
        # reorder tw_df using our ordered column names
        tw_df = tw_df[['label']+sorted_taxon_cols]
    else: # No labels. Just reorder by rank and name
        taxon_cols = tw_df.columns.tolist()
        rank_taxon = zip([rank_order.index(tw_df[c]['rank']) for c in taxon_cols],taxon_cols)
        sorted_taxon_cols = [rt[1] for rt in sorted(rank_taxon)]
        tw_df = tw_df[sorted_taxon_cols]
    
    # Flip so that taxons are rows, and columns are specimens
    tw_df_T = tw_df.T
    
    # Move rank to the left and sort by labels if they exist
    if args.labels:
        sorted_label_specimens = sorted([(tw_df_T[c]['label'],c) for c in tw_df_T.columns.tolist() if c != 'rank'])
        sorted_specimens = [ls[1] for ls in sorted_label_specimens]
        tw_df_T = tw_df_T[['rank']+sorted_specimens]
    else: # No labels
        sorted_specimens = sorted([c for c in tw_df_T.columns.tolist() if c != 'rank'])
        tw_df_T = tw_df_T[['rank']+sorted_specimens]
    
    
    # Output to CSV
    tw_df_T.to_csv(args.tallies_wide)
    
    # If we are also to make a fractions wide format output
    if args.fractions_wide:
        fw_df_T = pd.DataFrame()
        fw_df_T['rank'] = tw_df_T['rank']
        
        for specimen in sorted_specimens:
            if args.labels:
                fw_df_T[specimen] = tw_df_T[specimen].iloc[1:] / float(tw_df_T[specimen].iloc[1:].sum())
                fw_df_T[specimen]['label'] = tw_df_T.loc['label',specimen]
            else: # No labels, easier
                fw_df_T[specimen] = tw_df_T[specimen] / float(tw_df_T[specimen].sum())
        
        fw_df_T.to_csv(args.fractions_wide)
    

# Boilerplate method to run this as a script
if __name__ == '__main__':
    main()
