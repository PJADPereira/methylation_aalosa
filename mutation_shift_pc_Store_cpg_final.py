import os
import subprocess
import xlsxwriter
import multiprocessing
from itertools import combinations, permutations
from collections import defaultdict
import pickle

def parse_sfs(file_handle):
    dict_freqs = defaultdict(float)
    dict_bases = defaultdict(lambda:{"ref":None,"alt":None})


    for line in file_handle:
        try:
            chrm, pos, ref, nref,nalt,qref,qalt,bases,_,_,freq = line.strip().split()
        
        
        except ValueError:
            if line.strip().split()[3] == "*":
                continue
            else:
                print(line)
                input()

        else:
            
            if len(bases)==2:

                palt = bases.replace(ref,"")


                dict_freqs[f"{chrm}:{pos}"] = float(freq)
                
                if float(freq) > 0.5:

                    dict_bases[f'{chrm}:{pos}']['ref'] = palt
                    dict_bases[f'{chrm}:{pos}']['alt'] = ref

                else:
                    dict_bases[f'{chrm}:{pos}']['ref'] = ref
                    dict_bases[f'{chrm}:{pos}']['alt'] = palt


            elif len(bases)>2:
                print(line)
                input()

            else:
                continue

    return dict_freqs,dict_bases

def sanity_check(pos,p1b,p1f,p2b,p2f):
    if p1b[pos]['ref'] == p2b[pos]['alt']:
        if p2b[pos]['ref'] == p1b[pos]['alt']:
            return True
        else:

            # print(p1b[pos])
            # print(p2b[pos])
            # print(p1f[pos])
            # print(p2f[pos])
            return False
    else:
        # print(p1b[pos])
        # print(p2b[pos])
        # print(p1f[pos])
        # print(p2f[pos])
        return False

def cpg_type(rtri,dtri):
    # to do CpG types the reverse strand needs to be included
    # to do so we need to look at the first two bases of the triplet to make sure we are lookgitn at the proper place
    
    ms_rdi = rtri[1:]
    ms_ddi = dtri[1:]

    cs_rdi = rtri[0:2]
    cs_ddi = dtri[0:2]



    if ms_rdi == "CG" and ms_ddi == "TG":
        return 1
    
    elif cs_rdi == "CG" and cs_ddi == "CA":
        return 1

    elif ms_rdi == "TG" and ms_ddi =="CG":
        return 2
	
	## these now in rc
    elif cs_rdi == "CA" and cs_ddi =="CG":  
        return 2

    else:
        return 0


cache_samtools = defaultdict(str)
worksheets = []
folder = "./methyl_analysis/pcadapt/maf_aa/final_set/"
sfs = [x for x in os.listdir(folder) if x.endswith(".sfs")]

pairs = permutations(sfs,2)

results_book = xlsxwriter.Workbook('./methyl_analysis/mutation_shift/v1.3_corrected/mutations_results_Rc.xlsx')



for pop1,pop2 in pairs:
    comparison_name = f"{pop1[0:3]}-{pop2[0:3]}"
    print(comparison_name)
    with open(f"./methyl_analysis/mutation_shift/v1.3_corrected/{comparison_name}_cpgs_positions","w") as cpgs: 
        
        pop1 = f"{folder}{pop1}"
        pop2 = f"{folder}{pop2}"
        
        
        

        w1 = results_book.add_worksheet(f"{comparison_name}_1")
        w2 = results_book.add_worksheet(f"{comparison_name}_3")
        

        #load the files and parse both frequencies and bases
        with open(pop1) as pop1_handle:
            p1_freqs, p1_bases = parse_sfs(pop1_handle)

        with open(pop2) as pop2_handle:
            p2_freqs, p2_bases = parse_sfs(pop2_handle)



        shared_positions = set(p1_freqs).intersection(p2_freqs)


        mutations = defaultdict(int)
        mutations_tri = defaultdict(int)
        cpg_positions = defaultdict(str)
        errors = 0
        for pos in shared_positions:
            if abs(p1_freqs[pos] - p2_freqs[pos]) > 0.5:
                if sanity_check(pos,p1_bases,p1_freqs,p2_bases,p2_freqs):
                    mutations[f"{p1_bases[pos]['ref']}>{p1_bases[pos]['alt']}"] +=1
                    if pos not in cache_samtools.keys():

                        chrm,bpos = pos.split(":")

                        samtools = subprocess.run(['samtools','faidx','./GCF_017589495.1_AALO_Geno_1.1_genomic.fna',f"{chrm}:{int(bpos)-1}-{int(bpos)+1}"],capture_output=True)
                        samtools = samtools.stdout.decode().strip().split("\n")[-1]
                        cache_samtools[pos] = samtools.upper()

                    else:
                        samtools = cache_samtools[pos]

                    rtri =[x.upper() for x in samtools]
                    rtri[1] = p1_bases[pos]['ref']
                    rtri = "".join(rtri)

                    dtri =[x.upper() for x in samtools]
                    dtri[1] = p1_bases[pos]['alt']
                    dtri = "".join(dtri)

                    dtype = cpg_type(rtri,dtri)
                    
                    if dtype == 1:
                        cpg_positions[pos] = "CpG->TpG"
                        mutations["CpG->TpG"]+=1
                    elif dtype == 2:
                        cpg_positions[pos] = "TpG->CpG"
                        mutations["TpG->CpG"]+=1
                    mutations_tri[f"{rtri}>{dtri}"] +=1
                    
                
        row = 0
        col = 0
        sorted_keys = sorted(mutations.keys())
        for k in sorted_keys:
            w1.write(row,col,k)
            w1.write(row,col+1,mutations[k])
            row+=1
            
        row = 0
        col = 0
        sorted_keys = sorted(mutations_tri.keys())
        for k in sorted_keys:
            w2.write(row,col,k)
            w2.write(row,col+1,mutations_tri[k])
            row+=1


        sorted_keys = sorted(cpg_positions.keys(), key= lambda x: ( int(x.split(":")[0][-4:-2]), int(x.split(":")[1]) ))
        for k in sorted_keys:
            cpgs.write("{}\t{}\n".format("\t".join(k.split(":")),cpg_positions[k]))

results_book.close()

with open("./methyl_analysis/mutation_shift/v1.3_corrected/cache.pkl","wb") as pkl:
    pickle.dump(cache_samtools,pkl,protocol=pickle.HIGHEST_PROTOCOL)