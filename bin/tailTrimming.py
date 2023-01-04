#!/usr/bin/env python3
import argparse
import re

from Bio import SeqIO

def trim_reads(fileobj):

    pattern = (
        r'type=(?P<type>\w+),'
        r'start=(?P<start>\d+),'
        r'end=(?P<end>\d+),'
        r'length=(?P<length>\d+),'
        r'identity=(?P<identity>\d+(?:\.\d+)?)%'
    )

    for seq_record in SeqIO.parse(fileobj, "fasta"):
        match = re.search(pattern, seq_record.description)
        result = match.groupdict()

        start = int(result['start']) - 1
        end = int(result['end'])

        head = seq_record.seq[:start]
        tail = seq_record.seq[end:]

        seq_record.seq = head + tail

        yield seq_record

def get_argument_parser():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--fasta",
         type=str,
         required=True,
         help="input fasta file",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="output tail.trimmed.fasta file",
    )

    return parser

def main():

    parser = get_argument_parser()
    args = parser.parse_args()

    with open(args.fasta) as fasta, open(args.output, "w") as out:
        SeqIO.write(trim_reads(fasta), out, "fasta-2line")

if __name__ == "__main__":
    main()