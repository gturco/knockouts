#script to find perfect target region for gene
## look in exons 300-600 bps
## two or more exons interrupted by introns not as ideal
## blast 19 mer to check for off target hits
from flatfeature import Bed
from pyfasta import Fasta

def get_lens(exon):
     region = abs(exon[0] - exon[1])
     return region
     
def find_no_hits(pattern,gene_name,ref_fasta):
    no_matches = True
    remaining_keys = True
    print "pattern",pattern
    while no_matches and remaining_keys:
        for cds in ref_fasta.keys():
            if cds == gene_name: continue
            seq = ref_fasta[cds][:]
            match = seq.find(pattern)
            if match > 1:
               no_matches = False
        remaining_keys = False
    return no_matches

def blast_exon(gene,exon,f):
    gene_fasta = f[gene['accn']]
    region = range(exon[0] - gene['start'], exon[1] - gene['start'])
    unq_regions = list([(0,0)])
    consecutive = True
    print region, gene
    for start in region:
        stop = start + 19
        consecutive = unq_regions[-1][1] == (stop -1)
        if stop > region[-1] : continue
        pattern = gene_fasta[start:stop]
        no_hits = find_no_hits(pattern,gene['accn'],f)
        if no_hits and consecutive:
            new_end = (unq_regions[-1][0], stop)
            unq_regions = unq_regions[:-1]  + [new_end]
        elif no_hits:
            unq_regions.append((start,stop))
    region_sizes = map(get_lens,unq_regions)
    sort_regions = zip(region_sizes,unq_regions)
    sort_regions.sort()
    biggest_unq_region = sort_regions[-1]
    return biggest_unq_region

def main(gene_name,gene_bedfile, fastq_cds):
    b = Bed(gene_bedfile)
    f = Fasta(fastq_cds)
    gene = b.accn(gene_name)
    exons = gene["locs"]
    for exon in exons:
        exon_size = get_lens(exon)
        if exon_size >= 300:
            best_unq_region = blast_exon(gene,exon,f)
            #print best_unq_region






main("Os01g02110","ricetest.bed","ricetest.fasta")
            ### check length inbetween exons
            # len inbetween, ## hits, ##
            #(0)
            ### Choose exons with small lens inbetween
            ###
