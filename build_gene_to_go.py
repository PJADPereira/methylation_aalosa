from collections import defaultdict
import gzip

string_names = defaultdict(str)
gene_sizes = defaultdict(int)
with gzip.open("../accessorydata/STRG0020VDH.protein.info.v11.5.txt.gz","rt") as names:
    header = True
    for line in names:
        if header:
            header = False
            continue
        
        string,gene,size,description = line.strip().split("\t")
        string_names[string]=gene
        gene_sizes[string]=int(size)*3

terms = defaultdict(set)
with gzip.open("../accessorydata/STRG0020VDH.protein.enrichment.terms.v11.5.txt.gz","rt") as termsinput:
    c=0
    for line in termsinput:
        string,ftype,term,function = line.strip().split("\t")
        if ftype.endswith("(Gene Ontology)"):
            terms[string].add(term)


with open("gene_length.csv","w") as length:
    for k,v in gene_sizes.items():
        length.write("{}\t{}\n".format(string_names[k],v))

with open("goterms.csv","w") as termsout:
    for k,v in terms.items():
        for term in list(v):
            termsout.write("{}\t{}\n".format(string_names[k],term))
