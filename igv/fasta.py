import pysam


def get_data(fa_file, region=None):
    if None == region:
        with open(fa_file, "r") as f:
            return (f.read())
    else:
        if isinstance(region, str):
            region = parse_region(region)

        chr = region["chr"]
        start = region["start"] - 1
        end = region["end"]

        fasta = pysam.FastaFile(fa_file)
        slice_seq = fasta.fetch(chr, start, end)
        fasta.close()

        return slice_seq


class FastaReader:
    def __init__(self, path):
        self.fa_file = pysam.FastaFile(path)

    def slice(self, region=None):
        if None == region:
            with open(self.fa_file, "r") as f:
                return (f.read())
        else:
            if isinstance(region, str):
                region = parse_region(region)

            chr = region["chr"]
            start = region["start"] - 1
            end = region["end"]

            try:
                seq = self.fa_file.fetch(chr, start, end)
            except KeyError:
                chr = "chr" + chr
                seq = self.fa_file.fetch(chr, start, end)

            return seq


def parse_region(locus_str):
    tokens = locus_str.split(":")
    chr = tokens[0]

    t2 = tokens[1].split('-')
    start = int(t2[0].replace(',', ''))
    if len(t2) > 1:
        end = int(t2[1].replace(',', ''))
    else:
        end = start

    return {'chr': chr, 'start': start, 'end': end}
