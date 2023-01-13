# The following targets assume you have all of Baxolele's BAM files in */*.bam

all:
	for q in 13 20 30; \
        do \
            echo "Quality $$q"; \
            ./codon-count-with-soft-and-quality.py --minQuality $$q --site 27217 --bam samples/*/*.bam > out-qual-$$q.txt; \
        done

# This was my original processing, which did not include soft-clipped bases
# or exclude codons with any low-quality base.
original:
	for i in samples/*/*.bam; do echo "$$i"; ./codon-count.py --site 27217 "$$i"; echo; done
