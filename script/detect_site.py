import pandas as pd

site_aln_file = snakemake.input[0]
site_file = snakemake.output[0]
min_reads_num = snakemake.config['min_reads_num']
cluster_dist = snakemake.config['cluster_dist']

site_alns = pd.read_table(site_aln_file).to_dict('records')
aln_len = len(site_alns)

clusters = []
i = 0
while i < aln_len:
    cluster = [i]
    site_aln_i = site_alns[i]
    chr_h_i = site_aln_i['t_name_h']
    chr_v_i = site_aln_i['t_name_v']
    pos_h_i = site_aln_i['pos_h']
    pos_v_i = site_aln_i['pos_v']
    insert_len_i = site_aln_i['insert_len']

    j = i + 1
    while j < aln_len:
        site_aln_j = site_alns[j]
        chr_h_j = site_aln_j['t_name_h']
        chr_v_j = site_aln_j['t_name_v']
        pos_h_j = site_aln_j['pos_h']
        pos_v_j = site_aln_j['pos_v']
        insert_len_j = site_aln_j['insert_len']

        if chr_h_i == chr_h_j and chr_v_i == chr_v_j and \
            abs(pos_h_j - pos_h_i) <= cluster_dist and \
                abs(pos_v_j - pos_v_i) <= cluster_dist:
            cluster.append(j)
        else:          
            break

        j += 1

    if len(cluster) >= min_reads_num:
        clusters.append(cluster)
    
    i = j

fout = open(site_file, 'w')
fout.write('id\tchr_host\tposition_host\tchr_vector\tposition_vector\tread_num\trevmap\n')
for i in range(len(clusters)):
    cluster = clusters[i]
    pos_hosts = []
    pos_vectors = []
    site_aln = site_alns[cluster[0]]
    chr_host = site_aln['t_name_h']
    chr_vector = site_aln['t_name_v']
    
    for aln_id in cluster:
        site_aln = site_alns[aln_id]
        pos_hosts.append(site_aln['pos_h'])
        pos_vectors.append(site_aln['pos_v'])
    
    pos_hosts.sort()
    position_h = pos_hosts[0]

    pos_vectors.sort()
    position_v = pos_vectors[0]
    
    revmap = ','.join(str(i) for i in cluster)
    fout.write(f'{i}\t{chr_host}\t{position_h}\t{chr_vector}\t{position_v}\t{len(cluster)}\t{revmap}\n')

fout.close()