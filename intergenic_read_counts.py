import pysam
import gffutils
from collections import defaultdict

"""
Usage:
    This script identifies intergenic regions in a genome based on a provided GFF annotation file,
    limited to features that contain "ID=gene". It calculates the length of each intergenic region
    and determines the number of sequencing reads from a BAM alignment file that map to each region.

    The script performs the following steps:
    1. Parses the .gff file to identify all gene features and determines intergenic regions
       between adjacent genes on the same contig.
    2. Filters only features where the attribute string contains "ID=gene".
    3. Computes the length of each intergenic region.
    4. Uses the BAM file to count the number of aligned reads that overlap each intergenic region.

    Required Inputs:
        - A GFF annotation file (e.g., genome_annotations.gff)
        - A BAM alignment file (e.g., reads_aligned.bam)

    Example:
        python intergenic_read_counts.py genome_annotations.gff reads_aligned.bam

    Dependencies:
        - pysam
        - gffutils

    Output:
        Printed list of intergenic regions with coordinates, length, and read count:
        Format: Region: <seqid>:<start>-<end>, Length: <length>, Reads: <count>

"""

def parse_gff(gff_file):
    #print(f"Parsing GFF file: {gff_file}")
    # Create a GFF database from the file
    db = gffutils.create_db(gff_file, dbfn=':memory:', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
    #print("GFF database created in memory")

    # Extract intergenic regions based on gene features
    intergenic_regions = []
    genes = list(db.features_of_type('CDS'))
    #print(f"Number of gene features found: {len(genes)}")

    for i in range(len(genes) - 1):
        current_gene = genes[i]
        next_gene = genes[i + 1]
        
        if current_gene.seqid == next_gene.seqid:
            intergenic_start = current_gene.end + 1
            intergenic_end = next_gene.start - 1
            
            if intergenic_start < intergenic_end:
                region = {
                    'seqid': current_gene.seqid,
                    'start': intergenic_start,
                    'end': intergenic_end,
                    'length': intergenic_end - intergenic_start + 1,
                    'reads': 0
                }
                intergenic_regions.append(region)
                #print(f"Intergenic region added: {region['seqid']}:{region['start']}-{region['end']} (length {region['length']})")

    #print(f"Total intergenic regions identified: {len(intergenic_regions)}")
    return intergenic_regions

def count_reads_in_intergenic_regions(bam_file, intergenic_regions):
    #print(f"Opening BAM file: {bam_file}")
    bam = pysam.AlignmentFile(bam_file, "rb")

    for region in intergenic_regions:
        #print(f"Counting reads in region: {region['seqid']}:{region['start']}-{region['end']}")
        # Fetch reads that overlap the intergenic region
        reads = bam.fetch(region['seqid'], region['start'], region['end'])
        count = sum(1 for _ in reads)
        region['reads'] = count
        #print(f"Reads counted: {count}")

    bam.close()
    #print("BAM file closed")
    return intergenic_regions

def main(gff_file, bam_file):
    #print("Starting main process")
    intergenic_regions = parse_gff(gff_file)
    intergenic_regions = count_reads_in_intergenic_regions(bam_file, intergenic_regions)

    #print("Final intergenic region read counts:")
    for region in intergenic_regions:
        print(f"Region: {region['seqid']}:{region['start']}-{region['end']}, Length: {region['length']}, Reads: {region['reads']}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python script.py <annotations.gff> <alignments.bam>")
    else:
        gff_file = sys.argv[1]
        bam_file = sys.argv[2]
        main(gff_file, bam_file)
