# The following assumes you have all of Baxolele's BAM files in */*.bam
# It's just an example of how to call codon-count.py on a bunch of files.
all:
	for i in */*.bam; do echo "$$i"; ./codon-count.py --site 27217 "$$i"; echo; done
