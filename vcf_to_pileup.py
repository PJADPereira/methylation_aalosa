import sys
from collections import defaultdict
import pickle

vcf_path = sys.argv[1]#"./headed.vcf" sys.argv[1]
bed_path = sys.argv[2] #"final list of SNPs to consider" 


bed_snps = defaultdict(str)

with open(bed_path) as bed:
    has_header = True
    for line in bed:
        if has_header:
            has_header = False
            continue
        else:
            chrm,pos,anc_ref,*_ = line.strip().split()
            bed_snps[f"{chrm}:{pos}"] = anc_ref



quals = []
with open(vcf_path) as vcf:
    for line in vcf:
        header = line[1:].strip().split()
        break

chrm,pos,idd,ref,alt,qal,flter,info,formt,*pops = header
pops_order = pops
pops_snps = defaultdict(lambda:{x:list() for x in pops})



with open(vcf_path) as vcf:
    
    pileups = [open(f"{pop}_pileupfromvcf.pileup" ,"w") for pop in pops_order]
    has_header = True
    for line in vcf:

        if has_header:
            has_header = False
            continue
        else:
            chrm,pos,idd,ref,alt,qal,flter,info,formt,*pops = line.strip().split()
            bed_ref = bed_snps[f"{chrm}:{pos}"]
            if bed_ref:
                if bed_ref == "N":
                    continue

                elif len(set([bed_ref,ref,alt]))>2:
                    print("{}\t{}\ttriallelic".format(chrm,pos))   # UNCOMMENT ON THE FINAL VERSION
                    continue
                else:
                    for i, pop_data in enumerate(pops):
                        try:
                            _,dp,counts,_,rqal,_,aqal = pop_data.split(":")
                        except ValueError:
                            if bed_ref == ref:
                                continue
                                #pops_snps[f"{chrm}:{pos}"][pops_order[i]]=[0,0,0,0,ref,alt]
                            else:
                                continue
                                #pops_snps[f"{chrm}:{pos}"][pops_order[i]]=[0,0,0,0,alt,ref]

                        rqal = int(rqal)
                        aqal = int(aqal)
                        r,a = list(map(int,counts.split(",")))
                        #print(pops_order[i],chrm,pos)
                        perbaserqal = rqal//r if r >0 else 0
                        perbaseaqal = aqal//a if a >0 else 0
                        rqalstring = chr(33+perbaserqal)*r
                        aqalstring = chr(33+perbaseaqal)*a
                        if r+a == 0:
                            print("{}\t{}\tsnp was called based for this {} pop with just a single read that was ignored".format(chrm,pos,pops_order[i]]))   # UNCOMMENT ON THE FINAL VERSION
                            continue
                        if bed_ref == ref:
                            pileups[i].write("{}\t{}\t{}\t{}\t{}\t{}\n".format(chrm,pos,bed_ref,r+a,"."*r+alt*a,rqalstring+aqalstring))
                            #pops_snps[f"{chrm}:{pos}"][pops_order[i]]=[r,a,rqal,aqal,ref,alt]
                        elif bed_ref == alt:
                            pileups[i].write("{}\t{}\t{}\t{}\t{}\t{}\n".format(chrm,pos,bed_ref,r+a,"."*a+ref*r,aqalstring+rqalstring))
                            #pops_snps[f"{chrm}:{pos}"][pops_order[i]]=[a,r,aqal,rqal,alt,ref]
                        else:
                            print("{}\t{}\tweirdsnp that shoudlnt exist".format(chrm,pos))   # UNCOMMENT ON THE FINAL VERSION
                            continue
                        
            else:
                continue


[x.close() for x in pileups]
####pileup
#CHRM,pos,REF,DEPTH,bases,qual
#bases = . if equal to ref, else ALT allele

##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Number of observation for each allele">
##INFO=<ID=RO,Number=1,Type=Integer,Description="Count of full observations of the reference haplotype.">
##INFO=<ID=QR,Number=1,Type=Integer,Description="Reference allele quality sum in phred">
##INFO=<ID=AO,Number=A,Type=Integer,Description="Count of full observations of this alternate haplotype.">
##INFO=<ID=QA,Number=A,Type=Integer,Description="Alternate allele quality sum in phred">
