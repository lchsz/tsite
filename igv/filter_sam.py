import pandas as pd
import os
import util

site_file = snakemake.input[0]
site_aln_file = snakemake.input[1]
sam_v1_file = snakemake.input[2]
sam_v2_file = snakemake.input[3]
sam_h1_file = snakemake.input[4]
sam_h2_file = snakemake.input[5]
bam_dir = snakemake.output[0]

if not os.path.exists(bam_dir):
    os.mkdir(bam_dir)

sites = pd.read_table(site_file).to_dict('records')
site_alns = pd.read_table(site_aln_file, index_col ='id')

for site in sites:
    site_id = site['id']
    out_sam_file_v = os.path.join(bam_dir, str(site_id) + '_v.sam')
    out_bam_file_v = os.path.join(bam_dir, str(site_id) + '_v.bam')
    out_sam_file_h = os.path.join(bam_dir, str(site_id) + '_h.sam')
    out_bam_file_h = os.path.join(bam_dir, str(site_id) + '_h.bam')
    
    fout_v = open(out_sam_file_v, 'wt')
    fout_h = open(out_sam_file_h, 'wt')

    chr_names_v = set()
    chr_names_h = set()
    site_aln_ids = site['revmap'].split(',')
    for site_aln_id in site_aln_ids:
        site_aln = site_alns.loc[int(site_aln_id)]
        read_type = site_aln['read_type']
        chr_name_v = site_aln['t_name_v']
        chr_name_h = site_aln['t_name_h']

        if read_type == 'R1':
            if not chr_name_v in chr_names_v:
                util.filter_sam_header(chr_name_v, sam_v1_file, fout_v)
                chr_names_v.add(chr_name_v)
            if not chr_name_h in chr_names_h:
                util.filter_sam_header(chr_name_h, sam_h1_file, fout_h)
                chr_names_h.add(chr_name_h)
        elif read_type == 'R2':
            if not chr_name_v in chr_names_v:
                util.filter_sam_header(chr_name_v, sam_v2_file, fout_v)
                chr_names_v.add(chr_name_v)
            if not chr_name_h in chr_names_h:
                util.filter_sam_header(chr_name_h, sam_h2_file, fout_h)
                chr_names_h.add(chr_name_h)

    for site_aln_id in site_aln_ids:
        site_aln = site_alns.loc[int(site_aln_id)]
        read_name = site_aln['q_name']
        read_type = site_aln['read_type']

        if read_type == 'R1':
            util.filter_sam(read_name, sam_v1_file, fout_v)
            util.filter_sam(read_name, sam_h1_file, fout_h)
        elif read_type == 'R2':
            util.filter_sam(read_name, sam_v2_file, fout_v)
            util.filter_sam(read_name, sam_h2_file, fout_h)

    fout_v.close()
    fout_h.close()

    os.system(f"samtools view -b {out_sam_file_v} | samtools sort -o {out_bam_file_v} -")
    os.system(f"samtools index {out_bam_file_v}")
    os.system(f"samtools view -b {out_sam_file_h} | samtools sort -o {out_bam_file_h} -")
    os.system(f"samtools index {out_bam_file_h}")