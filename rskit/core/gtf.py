"""
Tools for working with GTF files.

GTF Format:

Fields must be tab-separated. Also, all but the final field in each feature line must contain a value; "empty" columns should be denoted with a '.'

    seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
    source - name of the program that generated this feature, or the data source (database or project name)
    feature - feature type name, e.g. Gene, Variation, Similarity
    start - Start position of the feature, with sequence numbering starting at 1.
    end - End position of the feature, with sequence numbering starting at 1.
    score - A floating point value.
    strand - defined as + (forward) or - (reverse).
    frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
    attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.

Some fields are renamed to match our common convention.

1   havana  gene    11869   14409   .   +   .   gene_id "ENSG00000223972"; gene_version "5"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene";

attr_format: setting the attribute format to 'ensembl' will use a simplified, and somewhat faster, parser for the attribute list, which is then expected to be well-formed.
             otherwise, the permissive regex described below is used.

"""

import re

re_attrs = re.compile(r'(\w+)(?:\s*=\s*|\s+)(?:"(.*?)"|(.*?));\s*')
id_keys = [
    ('gene_id', 'gene_version'),
    ('transcript_id', 'transcript_version'),
    ('protein_id', 'protein_version'),
]

# attribute tags are optionally quoted, and separated by wither whitespace or '='
#
# the match pattern is:
#  - word
#  - '=' OR whitespace
#  - '*' any '*' OR any (non-greedy match)
#  - ';'
#  - optional whitespace


def open(reader, attr_format=None, keep_line=False, append_versions=False):
    # strip only line terminators; by default, rstrip() will remove trailing tabs
    lines = (line.rstrip('\r\n') for line in reader if not line[0] == '#')
    for line in lines:
        rec = parse_line(line, attr_format, append_versions)
        yield (rec, line) if keep_line else rec


def parse_line(line, attr_format=None, append_versions=False):

    # 8 splits, 9 fields; attribute field may include tabs
    cols = [f if f != '.' else None for f in line.split('\t', 8)]

    rec = GtfRecord()
    rec.chr = cols[0]
    rec.source = cols[1]
    rec.feature = cols[2]
    rec.start = int(cols[3]) if cols[3] is not None else None
    rec.stop = int(cols[4]) if cols[4] is not None else None
    rec.score = float(cols[5]) if cols[5] is not None else None
    rec.strand = cols[6]
    rec.frame = int(cols[7]) if cols[7] is not None else None

    if attr_format == 'ensembl':
        rec.meta = dict((t, v[1:-1]) for t, v in (tag.split(' ', 1) for tag in cols[8].rstrip(';').split('; ')))

        if append_versions:
            keys = [k for k in rec.meta.keys() if k.endswith('_id')]
            for key, ver in id_keys:
                if key in rec.meta and ver in rec.meta:
                    rec.meta[key] = rec.meta[key] + '.' + rec.meta.pop(ver)

    else:
        rec.meta = dict((m[0], m[1] or m[2]) for m in re_attrs.findall(cols[8]))

    if 'gene_id' in rec.meta:
        rec.gene_id = rec.meta['gene_id']
        if rec.feature == 'gene':
            rec.id = rec.gene_id

    if 'transcript_id' in rec.meta:
        rec.transcript_id = rec.meta['transcript_id']
        if rec.feature == 'transcript':
            rec.id = rec.transcript_id

    return rec


class GtfRecord:

    def __init__(self):

        self.chr = None
        self.source = None
        self.feature = None
        self.start = None
        self.stop = None
        self.score = None
        self.strand = None
        self.frame = None
        self.meta = None

        self.id = None
        self.gene_id = None
        self.transcript_id = None

    def fields(self):
        fmeta = ' '.join(f'{k} "{v}";' for k, v in self.meta.items())
        return [self.chr, self.source, self.feature, self.start, self.stop, self.score, self.strand, self.frame, fmeta]

    def __str__(self):
        return '\t'.join(str(f if f is not None else '.') for f in self.fields())
