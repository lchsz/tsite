import pysam


def build_aliastable(chrs):
    chralias = {}

    for c in chrs:
        if c.startswith("chr"):
            if c == 'chrM':
                chralias['MT'] = 'chrM'
            else :
                alias = c[3:]
                chralias[alias] = c
        else:
            if c == 'MT':
                chralias['chrM'] = 'MT'
            else:
                alias = f'chr{c}'
                chralias[alias] = c

    return chralias


def get_alias(c):
    if c.startswith("chr"):
        if c == 'chrM':
            return 'chrM'
        else :
            return c[3:]
    else:
        if c == 'MT':
            return 'MT'
        else:
            return f'chr{c}'


class BamReader:
    def __init__(self, filename, args=None):
        self.filename = filename
        self.args = args

        samargs = ["-H", filename]
        header = pysam.view(*samargs)
        seqnames = parse_seqnames(header)
        self.aliastable = build_aliastable(seqnames)

    # add sam flag for unit tests
    def slice(self, region=None, region2=None,  sam=False):
        if sam:
            samargs = ["-h", self.filename]
        else:
            samargs = ["-b", "-h", self.filename]

        samargs.append("-F")
        samargs.append("1536")

        if region:
            range_string = self.get_chrname(
                region['chr']) + ":" + str(region['start']) + "-" + str(region['end'])
            samargs.append(range_string)
        if region2:
            range_string = self.get_chrname(
                region2['chr']) + ":" + str(region2['start']) + "-" + str(region2['end'])
            samargs.append(range_string)

        return pysam.view(*samargs)

    def get_chrname(self, c):
        if c in self.aliastable:
            return self.aliastable[c]
        else:
            return c
        

def parse_seqnames(header):
    seqnames = []
    lines = header.split('\n')
    for line in lines:
        if line.startswith('@SQ'):
            idx1 = line.find("SN:")
            if idx1 > 0:
                idx2 = line.find("\t", idx1)
                seqnames.append(line[idx1+3:idx2])

    return seqnames
