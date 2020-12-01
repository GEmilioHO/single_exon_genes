#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 12:39:33 2020

@author: andresgarcíagarcía
@edited: gabrielemilioherreraoropeza

"""

import sys
import ensembl_rest
import fetch_genes


def single_exon_isoforms(gene_info):
    return ','.join([transcript['id'] for transcript in gene_info['Transcript']
                             if len(transcript['Exon']) == 1])
# ---

species = (sys.argv[1]
            if (len(sys.argv) > 1) 
            else 'gallus gallus')

print("Species:", species)

client = ensembl_rest.EnsemblClient()

# --- Fetch information of the chromosomes in Gallus gallus
print('Fetching chromosome data...')

genome_info = client.assembly_info(species=species)
chrom_data = {item['name']: item['length'] 
                for item in genome_info['top_level_region']
                if item['coord_system'] == 'chromosome'} 

# --- Fetch the ids of all the genes in the Gallus Gallus genome
print('Fetching all gene IDs...')

all_genes = sum((fetch_genes.genes_in_chrom(species, chrom, length) 
                    for chrom, length in chrom_data.items()), 
                []
            )


# --- Get the list of exons for all genes
print('Fetching gene info...')

fields = ['id', 'display_name', 'description', 'biotype', 'start', 'end', 'seq_region_name']

genes = fetch_genes.get_info(all_genes, 
                             (lambda gene_info : {
                                 "exons_n": len({exon['id'] 
                                                 for transcript in gene_info['Transcript']
                                                 for exon in transcript['Exon']}),
                                 "transcripts_n": len(gene_info['Transcript']),
                                 "single_exon_isoforms": single_exon_isoforms(gene_info),
                                 **{field:gene_info[field] for field in fields if field in gene_info}
                                 }
                             )
        )

# --- Filter those which only have one exon ---
single_exon_genes = {gene_id for gene_id, gene_info in genes.items() if gene_info['exons_n'] == 1}
single_exon_isoforms = {gene_id for gene_id, gene_info in genes.items() 
                            if gene_info['single_exon_isoforms'] and gene_id not in single_exon_genes}
coding_genes = {gene_id for gene_id, gene_info in genes.items() if gene_info['biotype'] == 'protein_coding'}


# --- Report the results
print(f"Of the {len(all_genes)} genes in {species}, "
      f"only {len(single_exon_genes)} are single exon genes, "
      f"{len(coding_genes)} are protein coding,"
      f"and {len(single_exon_isoforms)} multi exon genes have single exon isoforms.")

# --- Write the list of the single exon gene IDs
import csv

with open(f"../data/{species}-single_exon_genes.tsv", 'w', newline='') as seg_out, \
     open(f"../data/{species}-single_exon_isoforms.tsv", 'w', newline='') as sei_out, \
     open(f"../data/{species}-single_exon_genes_protein_coding.tsv", 'w', newline='') as seg_cod_out:

        fieldnames = fields + ['exons_n', 'transcripts_n', 'single_exon_isoforms']
        seg_writer = csv.DictWriter(seg_out, fieldnames=fieldnames, dialect='excel-tab')
        sei_writer = csv.DictWriter(sei_out, fieldnames=fieldnames, dialect='excel-tab')
        seg_cod_writer = csv.DictWriter(seg_cod_out, fieldnames=fieldnames, dialect='excel-tab')

        seg_cod_writer.writeheader()
        seg_writer.writeheader()
        sei_writer.writeheader()
        for gene_id, gene_info in genes.items():
            if gene_id in coding_genes and gene_id in single_exon_genes:
                seg_cod_writer.writerow(gene_info)            
            if gene_id in single_exon_genes:
                seg_writer.writerow(gene_info)
            elif gene_id in single_exon_isoforms:
                sei_writer.writerow(gene_info)