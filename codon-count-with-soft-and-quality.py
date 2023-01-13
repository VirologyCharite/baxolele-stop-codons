#!/usr/bin/env python

import sys
import argparse
from collections import defaultdict
from operator import itemgetter

from dark.aa import STOP_CODONS
from dark.sam import samfile
from dark.utils import pct


DEFAULT_REF_ID = "hCoV-19/Wuhan-Hu-1/2019|EPI_ISL_402125"


def makeParser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Find stop codons for Baxolele")

    parser.add_argument(
        "--bam",
        nargs="+",
        metavar="FILE.bam",
        help="The BAM file(s) to read.",
    )

    parser.add_argument(
        "--site",
        metavar="N",
        required=True,
        type=int,
        help=(
            "The 1-based site location of the first nucleotide in the codon "
            "to check."
        ),
    )

    parser.add_argument(
        "--minQuality",
        metavar="N",
        type=int,
        default=30,
        help="The minimum acceptable quality (PHRED) score.",
    )

    parser.add_argument(
        "--reference",
        metavar="SEQ-NAME",
        help="The name of sequence in the BAM file to find matches for.",
    )

    return parser


def adjustInitialReferenceOffsets(offsets):
    result = offsets[:]
    leadingNoneCount = 0
    firstNonNoneOffset = None

    for offset in result:
        if offset is None:
            leadingNoneCount += 1
        else:
            firstNonNoneOffset = offset
            break

    for i in range(leadingNoneCount):
        result[i] = -1 * (firstNonNoneOffset - leadingNoneCount + i)

    return result


# Test
assert adjustInitialReferenceOffsets([None, None, 3, 4]) == [-1, -2, 3, 4]
assert adjustInitialReferenceOffsets([None, None, None, 5, 6, 7]) == [
    -2,
    -3,
    -4,
    5,
    6,
    7,
]


def adjustFinalReferenceOffsets(offsets):
    result = offsets[::-1]
    trailingNoneCount = 0
    lastNonNoneOffset = None

    for offset in result:
        if offset is None:
            trailingNoneCount += 1
        else:
            lastNonNoneOffset = offset
            break

    for i in range(trailingNoneCount):
        result[i] = -1 * (lastNonNoneOffset + trailingNoneCount - i)

    return result[::-1]


# Test
assert adjustFinalReferenceOffsets([3, 4, None, None]) == [3, 4, -5, -6]
assert adjustFinalReferenceOffsets([5, 6, 7, None, None, None]) == [
    5,
    6,
    7,
    -8,
    -9,
    -10,
]


def processFile(
    filename: str, reference: str, offset: int, minQuality: int, wantedOffsets: set[int]
) -> None:
    with samfile(filename) as bam:
        stopCount = readCount = lowQualityReadCount = nonOverlappingReadCount = 0
        partiallyOverlappingCount = 0
        codons: dict[str, int] = defaultdict(int)
        ids = set()
        lowQualityIds = set()
        nonOverlappingIds = set()
        partiallyOverlappingIds = set()
        for readCount, read in enumerate(
            bam.fetch(reference, start=offset, stop=offset + 3),
            start=1,
        ):
            id_ = read.query_name
            if read.is_read1:
                id_ += "-1"
            else:
                id_ += "-2"
                assert read.is_read2
            assert id_ not in ids
            ids.add(id_)
            assert not read.is_unmapped
            bases = [None, None, None]
            foundCount = 0

            referenceOffsets = adjustFinalReferenceOffsets(
                adjustInitialReferenceOffsets(
                    list(read.get_reference_positions(full_length=True))
                )
            )

            # Make sure the reference offsets include some of the wanted offsets.
            if not set(referenceOffsets) & wantedOffsets:
                # print("Reference offsets", read.cigarstring, read.get_blocks(),
                #       "do not overlap wanted offsets", file=sys.stderr)
                nonOverlappingReadCount += 1
                nonOverlappingIds.add(id_)
                continue

            # Make sure there are no duplicates in the reference offsets.
            # This is too simplistic, since None can be in the list multiple times,
            # so commented out for now.
            # assert len(set(referenceOffsets)) == len(referenceOffsets)

            querySequence = list(read.query_sequence)
            queryQualities = list(read.query_qualities)

            assert len(querySequence) == len(queryQualities) == len(referenceOffsets)

            for queryBase, queryQuality, referenceOffset in zip(
                querySequence, queryQualities, referenceOffsets
            ):
                if referenceOffset is None:
                    # This is a reference offset that is missed by the read
                    # (i.e., a deletion in the read).
                    continue

                if referenceOffset < 0:
                    referenceOffset *= -1
                    soft = True
                else:
                    soft = False

                if referenceOffset in wantedOffsets:
                    if queryQuality < minQuality:
                        # There is a base in the target region with a
                        # too-low quality. Skip the read.
                        assert foundCount < 3
                        lowQualityReadCount += 1
                        lowQualityIds.add(id_)
                        break

                    codonOffset = referenceOffset - offset
                    assert bases[codonOffset] is None
                    bases[codonOffset] = queryBase

                    foundCount += 1
                    if foundCount == 3:
                        # All codon offsets have been seen, so we can stop
                        # processing this read.
                        break

            if foundCount == 3:
                codon = "".join(bases)
                codons[codon] += 1
                stopCount += codon in STOP_CODONS
            else:
                partiallyOverlappingCount += 1
                partiallyOverlappingIds.add(id_)

        # assert not (partiallyOverlappingIds & lowQualityIds)
        assert not (partiallyOverlappingIds & nonOverlappingIds)
        assert not (lowQualityIds & nonOverlappingIds)

        print(f"{bam.mapped:,} reads were mapped to the reference.")
        print(
            f"{pct(readCount, bam.mapped)} overlap the region "
            f"({offset + 1}-{offset + 3}) in some way."
        )
        print("Of these:")
        print(
            f"  {pct(lowQualityReadCount, readCount)} have a low-quality "
            f"(<{minQuality}) base in the target region (ignored)."
        )
        print(
            f"  {pct(partiallyOverlappingCount, readCount)} only partially "
            f"overlap the region (ignored)."
        )
        print(
            f"  {pct(nonOverlappingReadCount, readCount)} span the region "
            f"but with no nucleotides matching in the specific sites (ignored)."
        )

        readsWithCodon = readCount - len(
            lowQualityIds | nonOverlappingIds | partiallyOverlappingIds
        )

        print(f"  {pct(stopCount, readsWithCodon)} of reads with a valid codon have a stop.")

        print(f"  Frequency breakdown of {readsWithCodon} reads with a valid codon:")

        total = 0
        for item in sorted(codons.items(), reverse=True, key=itemgetter(1)):
            codon, count = item
            total += count
            suffix = "\t*" if codon in STOP_CODONS else ""
            print(f"    {codon} {count:5d}{suffix}")

        if total != readsWithCodon:
            print(f"Huh!? {total} != {readsWithCodon} in {filename!r}", file=sys.stderr)
            sys.exit(1)


def main() -> None:
    args = makeParser().parse_args()
    offset = args.site - 1
    minQuality = args.minQuality
    wantedOffsets = set(range(offset, offset + 3))
    reference = args.reference or DEFAULT_REF_ID
    manyFiles = len(args.bam) > 1

    for filename in args.bam:
        if manyFiles:
            print(f"Processing {filename!r}", file=sys.stderr)
            print(f"Results for file {filename!r}")
        processFile(filename, reference, offset, minQuality, wantedOffsets)
        if manyFiles:
            print()


if __name__ == "__main__":
    main()
