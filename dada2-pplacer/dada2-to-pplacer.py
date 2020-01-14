#!/usr/bin/env python
import pandas as pd
import csv
import argparse
import sys

def main():
    args_parser = argparse.ArgumentParser()
    
    args_parser.add_argument('--seqtable','-s', help="Sequence table from dada2, in CSV format", required =True, type=file)
    args_parser.add_argument('--fasta_out_sequences', '-f', help="Write sequence variants to this file, in FASTA format", type=argparse.FileType('w'))
    args_parser.add_argument('--map', '-m', help="Write pplacer-style mapping of otu to specimen", type=argparse.FileType('w'))
    args_parser.add_argument('--weights', '-w', help="Write pplacer-style weights of otu by specimen", type=argparse.FileType('w'))
    
    args = args_parser.parse_args()
    
    # Check to see if we've been tasked with anything. If not, we have nothing to do and should exit
    if not (args.fasta_out_sequences or args.map or args.weights):
        sys.exit("Nothing to do")
    
    # Just convert our handles over to soething nicer
    
    if args.fasta_out_sequences:
        out_otu_seqs_h = args.fasta_out_sequences
    else:
        out_otu_seqs_h = None
    if args.map:
        out_map_h   = args.map
        map_writer = csv.writer(out_map_h)
    else:
        map_writer = None
    if args.weights:
        out_weights_h = args.weights
        weights_writer = csv.writer(out_weights_h)
    else:
        weights_writer=None
    
    # Load the sequence table
    seqtab = pd.read_csv(args.seqtable, index_col=0)
    # Generate OTU labels for each sequence variant, and create a dictionary that maps sequence-variant to otu_id...
    seq_to_otu_num = {s: 'otu-%d' % i for i, s in enumerate(seqtab.columns) }
    # and a dictionary to map otu_id to sequence-variant
    otu_num_to_seq = {'otu-%d' % i: s for i, s in enumerate(seqtab.columns) }
    # rename our columns to fit these new mappings
    seqtab.columns = [seq_to_otu_num[s] for s in seqtab.columns]
    
    # Annoyingly, we need to pick a representitive actual sequence from each OTU to be it's champion for guppy. To do so, we will go through each column, find the max count for that OTU, and use that specimen as the champion
    
    max_spec_for_otu = {otu_id: spec for otu_id, spec in seqtab.apply(lambda c: c.idxmax()).iteritems()}

    if out_otu_seqs_h:
        # Write out the sequences in fasta format, using the otu-id's generated above as an ID
        for otu_id, seq in otu_num_to_seq.iteritems():
            out_otu_seqs_h.write(">%s:%s\n%s\n" % (otu_id, max_spec_for_otu[otu_id] ,seq))

    # Now write the mapping and weights files
    # Both are headerless CSV format files
    # map: sequence_id (otu_id:specimen), specimen
    # weight: sequence_id (otu_id here), specimen_sequence_id (otu_id:specimen here), count
    # This is a bit of a clunky structure (relating to some historic cruft)

    if map_writer or weights_writer:
        for spec, row in seqtab.iterrows():
            row_nonzero = row[row>0]
            for otu_id, count in row_nonzero.iteritems():
                if map_writer:
                    map_writer.writerow([str(otu_id)+":"+str(spec),spec])
                if weights_writer:
                    weights_writer.writerow([otu_id+":"+max_spec_for_otu[otu_id],str(otu_id)+":"+str(spec),count])

    # Cleanup.    
    if out_otu_seqs_h:
        out_otu_seqs_h.close()
    if map_writer:
        out_map_h.close()
    if weights_writer:
        out_weights_h.close()

# Boilerplate method to run this as a script
if __name__ == '__main__':
    main()
