
from itertools import islice
import ensembl_rest
import sys
from tqdm import tqdm


def head(iterable, n=10):
    "Get the first entries from an iterable."
    iterator = iter(iterable)
    yield from islice(iterator, n)
# ---
        
def chunks_of(iterable, chunk_size=10):
    "Get the entries of the iterable in chunks."
    iterator = iter(iterable)
    while True:
        next_ = list(head(iterator, chunk_size))
        if next_:
            yield next_
        else:
            break
# ---



_client = ensembl_rest.EnsemblClient()

def unprocessed(arg): return arg

def overlapping_features(species, region, 
                         feature='gene', 
                         process=unprocessed, 
                         ensembl_client=_client):
    """Get the ids of all genes overlapping the given region.
    """
    raw_data = ensembl_client.overlap_region(species=species, 
                                             region=region, 
                                             params={'feature':feature})
    
    return [process(feature_data) for feature_data in raw_data]
# ---

def divide_region(start, end, max_region_length):
    # The Ensembl API only allows to query regions of 
    # length at most 5Mb at a time
    subregions = []
    region_upper = start
    while region_upper < end:
        new_region = (region_upper+1, 
                      region_upper + min(max_region_length, end-region_upper))
        subregions.append(new_region)
        
        region_upper += max_region_length

    return subregions
# ---

def genes_in_chrom(species, 
                   chrom, 
                   chrom_length=None, 
                   ensembl_client=_client,
                   max_region_length=5_000_000):
    """
    Fetch info of all the genes in the chromosome.
    """
    if chrom_length is None:
        chrom_length = ensembl_client.assembly_stats(species=species, 
                                                     region=chrom    )['length']
    
    # The Ensembl API only allows to query regions of 
    # length at most 5Mb at a time
    regions_to_query = divide_region(0, chrom_length, max_region_length)
        
    genes = []
    for start, end in regions_to_query:
        region = ensembl_rest.region_str(chrom, start, end)

        genes.extend(overlapping_features(species, region, 
                                          process=(lambda gene_info: gene_info['id']),
                                          ensembl_client=ensembl_client))
    
    return genes
# ---


def get_info(feature_ids, process_data, features_per_query=100, ensembl_client=_client):
    """Get the list of exons for each gene from a list of gene ids.
    """
    info = {}
    processed_features = 0
    for chunk in tqdm(chunks_of(feature_ids, features_per_query)):

        #print(processed_features)

        new_data = ensembl_client.lookup_post(params={'ids': chunk, 
                                                      'expand': True})

        info.update({gene_id:process_data(gene_data)
            for gene_id,gene_data in new_data.items()
            if gene_data
        })

        processed_features += len(chunk)
        
        
    return info
# ---