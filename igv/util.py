import os
from gzip import compress
from base64 import b64encode
import math
import json
import bam


def get_data_uri(data):
    '''
    Return a data uri for the input, which can be either a string or byte array
    '''
    if isinstance(data, str):
        data = compress(data.encode())
        mediatype = 'data:application/gzip'
    else:
        if data[0] == 0x1f and data[1] == 0x8b:
            mediatype = 'data:application/gzip'
        else:
            mediatype = 'data:application:octet-stream'

    enc_str = b64encode(data)
    data_uri = mediatype + ';base64,' + str(enc_str)[2:-1]
    return data_uri


def get_track_json(track_file):
    name = os.path.basename(track_file)
    return {
        'name': name,
        'url': track_file,
        'type': 'alignment',
        'format': 'bam',
        'showSoftClips': 'true'
    }


def locus_str(chr, start, end):
    if (end - start) == 1:
        return f'{chr}:{start + 1}'
    else:
        return f'{chr}:{start + 1}-{end}'


def load_session_uri(chr, position, flanking, fasta_reader, track):
    session = {}
    start = int(math.floor(position - flanking / 2))
    start = max(start, 1)
    end = int(math.ceil(position + flanking / 2))
    region = {
        'chr': chr,
        'start': start,
        'end': end
    }
    # Fasta
    seq = fasta_reader.slice(region)
    fa = '>' + chr + ':' + str(start) + '-' + str(end) + '\n' + seq
    fasta_uri = get_data_uri(fa)
    fasta_json = {
        'fastaURL': fasta_uri,
    }
    # Initial locus
    initial_locus = locus_str(chr, start, end)
    session = {
        'locus': initial_locus,
        'reference': fasta_json,
        'tracks': []
    }

    bam_reader = bam.BamReader(track['url'])
    bam_data = bam_reader.slice(region)
    track['url'] = get_data_uri(bam_data)
    session['tracks'].append(track)
    session_uri = get_data_uri(json.dumps(session))

    return session_uri


def filter_sam_header(chr_name, in_file, fout):
    fin = open(in_file, 'rt')
    for line in fin:
        line = line.strip()
        if len(line) == 0:
            continue

        ss = line.split('\t')
        if line.startswith('@SQ'):
            if ss[1] == 'SN:' + chr_name:
                fout.write(line + '\n')
    fin.close()


def filter_sam(read_name, in_file, fout):
    fin = open(in_file, 'rt')
    for line in fin:
        line = line.strip()
        if len(line) == 0:
            continue

        ss = line.split('\t')
        if read_name == ss[0]:
            fout.write(line + '\n')

    fin.close()
