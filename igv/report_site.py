from fasta import FastaReader
import json
from jinja2 import Environment, select_autoescape, FileSystemLoader
import pandas as pd
import os
import util

site_file = snakemake.input['site_file']
site_aln_file = snakemake.input['site_aln_file']
host_file = snakemake.input['host_file']
vector_file = snakemake.input['vector_file']
bam_dir = snakemake.input['bam_dir']
output_file = snakemake.output[0]

flanking = 300

sites = pd.read_table(site_file).to_dict('records')
site_alns = pd.read_table(site_aln_file).to_dict('records')
fa_reader_h = FastaReader(host_file)
fa_reader_v = FastaReader(vector_file)

sessions_h = {}
sessions_v = {}
id2revmap = {}
id2site_aln = {}

for site_aln in site_alns:
    id2site_aln[site_aln['id']] = site_aln

for site in sites:
    id2revmap[site['id']] = site['revmap'].split(',')
    session_id = str(site['id'])
    chr_h = site['chr_host']
    position_h = site['position_host']
    chr_v = site['chr_vector']
    position_v = site['position_vector']
    bam_file_h = os.path.join(bam_dir, f'{session_id}_h.bam')
    track_h = util.get_track_json(bam_file_h)
    session_uri_h = util.load_session_uri(chr_h, position_h, flanking,
                                          fa_reader_h, track_h)
    sessions_h[session_id] = session_uri_h
    bam_file_v = os.path.join(bam_dir, f'{session_id}_v.bam')
    track_v = util.get_track_json(bam_file_v)
    session_uri_v = util.load_session_uri(chr_v, position_v, flanking,
                                          fa_reader_v, track_v)
    sessions_v[session_id] = session_uri_v

sessions_v = json.dumps(sessions_v)
sessions_h = json.dumps(sessions_h)
id2revmap = json.dumps(id2revmap)
id2site_aln = json.dumps(id2site_aln)

loader = FileSystemLoader(searchpath='igv/templates')
env = Environment(
    loader=loader,
    autoescape=select_autoescape()
)
tpl = env.get_template('template_site.jinja2')
output = tpl.render(sites=sites, id2revmap=id2revmap, id2site_aln=id2site_aln,
                    sessions_v=sessions_v, sessions_h=sessions_h)
with open(output_file, 'w') as fout:
    fout.write(output)
