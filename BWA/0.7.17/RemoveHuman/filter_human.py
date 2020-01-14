from Bio import SeqIO
import argparse
import sys
import pysam
import gzip

###
#   For shotgun microbiome experiments, we wish to strip out human
#   genomic DNA (for various reasons, including privacy concerns).
#   A typical pipeline will align a human reference genome to remove any
#   reads from human genomic DNA sequences.
#   
#   This utility will take a SAM alignment file from paired end reads 
#   and filter the original read FASTQ files do those reads without
#   high-likelihood alignments to human
#
###


def get_passing_ids(aln_h, CUTOFFMAPQ, COV_MIN):
    r_1_aln = None
    passed_ids = set()
    with pysam.AlignmentFile(aln_h, 'r') as sf:
        for aln in sf:
            if aln.is_read1 is True:
                r_1_aln = aln
                continue
            # Implicit else we are read 2
            # First check to be sure we have the proper pair cached.
            assert(aln.qname == r_1_aln.qname)
            if (aln.is_unmapped or r_1_aln.is_unmapped):
                passed_ids.add(aln.qname)
            elif r_1_aln.alen / float(r_1_aln.query_length) <= COV_MIN or aln.alen / float(aln.query_length) <= COV_MIN:
                passed_ids.add(aln.qname)
            elif aln.mapq <= CUTOFFMAPQ or r_1_aln.mapq < CUTOFFMAPQ:
                passed_ids.add(aln.qname)
            else:  # Both forward and reverse reads seem to be a hit for human.
                continue
    return passed_ids


def filter_fastq(
        passed_ids,
        R1_fn,
        R2_fn,
        R1_out_fn,
        R2_out_fn):
    # Input files. Use FN to see if we need gzip
    if R1_fn[-2:].lower() == 'gz':
        R1_h = gzip.open(R1_fn, 'rt')
    else:
        R1_h = open(R1_fn, 'rt')
    if R2_fn[-2:].lower() == 'gz':
        R2_h = gzip.open(R2_fn, 'rt')
    else:
        R2_h = open(R2_fn, 'rt')

    # Output files. Use FN to see if we need gzip
    if R1_out_fn[-2:].lower() == 'gz':
        R1_out_h = gzip.open(R1_fn, 'wt')
    else:
        R1_out_h = open(R1_out_fn, 'wt')
    if R2_out_fn[-2:].lower() == 'gz':
        R2_out_h = gzip.open(R2_out_fn, 'wt')
    else:
        R2_out_h = open(R2_out_fn, 'wt')

    for R1, R2 in zip(
            SeqIO.parse(R1_h, 'fastq'),
            SeqIO.parse(R2_h, 'fastq')):
        assert(R1.id == R2.id)
        if R1.id in passed_ids:
            SeqIO.write(R1, R1_out_h, 'fastq')
            SeqIO.write(R2, R2_out_h, 'fastq')


def main():
    """Entrypoint for main script."""
    arg_parser = argparse.ArgumentParser(description="""
    This utility will take a SAM alignment file from paired end reads 
    and filter the original read FASTQ files do those reads without
    high-likelihood alignments to human.
    For gzipped alignments, consider using pipes: 
    gunzip -c ref.fna.gz | strip_mt_ebv.py | gzip > ref.nomtebv.fna.gz
    """)

    arg_parser.add_argument(
        '--alnfile', '-A',
        type=argparse.FileType('r'),
        help='Alignment File. Can be stdin. For gzip, consider pipes',
        default=sys.stdin
    )
    arg_parser.add_argument(
        '--r1in', '-1',
        required=True,
        help='Input fastq file for R1'
    )
    arg_parser.add_argument(
        '--r2in', '-2',
        required=True,
        help='Input fastq file for R2'
    )
    arg_parser.add_argument(
        '--r1out', '-o1',
        required=True,
        help='Output fastq file for R1'
    )
    arg_parser.add_argument(
        '--r2out', '-o2',
        required=True,
        help='Output fastq file for R2'
    )
    arg_parser.add_argument(
        '--mapq',
        default=30,
        type=int,
        help='Minimum mapq required to be considered a valid read'
    )
    arg_parser.add_argument(
        '--cov_min',
        type=float,
        default=0.9
    )

    args = arg_parser.parse_args()

    passed_ids = get_passing_ids(
        args.alnfile,
        args.mapq,
        args.cov_min,
    )

    filter_fastq(
        passed_ids,
        args.r1in,
        args.r2in,
        args.r1out,
        args.r2out
    )

if __name__ == "__main__":
    main()
