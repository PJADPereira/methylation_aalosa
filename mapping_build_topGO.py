from collections import defaultdict

mappings = defaultdict(list)

with open("./STRINGS/src/goterms.csv") as f:
    for line in f:
        name , term = line.strip().split()
        mappings[name].append(term)


with open("./categories_topGO.csv","w") as f:
    for k,v in mappings.items():
        f.write("{}\t{}\n".format(k,",".join(v)))