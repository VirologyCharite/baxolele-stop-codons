#!/usr/bin/env python

import argparse
from collections import defaultdict

from dark.aa import STOP_CODONS
from dark.sam import samfile
from dark.utils import pct


DEFAULT_REF_ID = "hCoV-19/Wuhan-Hu-1/2019|EPI_ISL_402125"


def makeParser() -> argparse.ArgumentParser:
    _parser = argparse.ArgumentParser(description="Find stop codons for Baxolele")

    _parser.add_argument("bam", metavar="FILE.bam", help="The BAM file to read.")

    _parser.add_argument(
        "--site",
        metavar="N",
        required=True,
        type=int,
        help=(
            "The 1-based site location of the first nucleotide in the codon "
            "to check."
        ),
    )

    _parser.add_argument(
        "--reference",
        metavar="SEQ-NAME",
        help="The name of sequence in the BAM file to find matches for.",
    )

    return _parser


def main() -> None:
    args = makeParser().parse_args()
    offset = args.site - 1
    wantedOffsets = set(range(offset, offset + 3))
    stopCount = 0
    codons = defaultdict(int)

    with samfile(args.bam) as bam:
        readCount = 0
        for readCount, read in enumerate(
            bam.fetch(args.reference or DEFAULT_REF_ID, offset, offset + 3), start=1
        ):
            assert not read.is_unmapped
            seen = set()
            query = read.query_sequence
            bases = [None, None, None]
            foundCount = 0
            for queryOffset, referenceOffset in read.get_aligned_pairs():
                if queryOffset is None or referenceOffset is None:
                    continue
                assert referenceOffset not in seen
                seen.add(referenceOffset)
                if referenceOffset in wantedOffsets:
                    codonOffset = referenceOffset - offset
                    assert bases[codonOffset] is None
                    bases[codonOffset] = query[queryOffset]
                    foundCount += 1
                    if foundCount == 3:
                        break

            if foundCount == 3:
                codon = "".join(bases)
                codons[codon] += 1
                stopCount += codon in STOP_CODONS

        print(f"Found {pct(stopCount, readCount)} stop codon reads")

        for codon, count in sorted(codons.items()):
            print(codon, count)


if __name__ == "__main__":
    main()
