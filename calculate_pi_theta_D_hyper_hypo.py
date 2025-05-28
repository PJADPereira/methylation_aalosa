import os
import subprocess

pool_sizes = {
"MON_AA_pileupfromvcf.pileup": 64,
"BAG_AA_pileupfromvcf.pileup": 58,
"CBD_AA_pileupfromvcf.pileup": 42,
"ALQ_AA_pileupfromvcf.pileup": 132,
"GAR_AA_pileupfromvcf.pileup": 96
}


max_coverage = {
"MON_AA_pileupfromvcf.pileup" : 74,
"BAG_AA_pileupfromvcf.pileup" : 72,
"CBD_AA_pileupfromvcf.pileup" : 54,
"ALQ_AA_pileupfromvcf.pileup" : 78,
"GAR_AA_pileupfromvcf.pileup" :60
}


fpool_sizes = {
"MONAASAB3_filtered.pileup": 64,
"BAGAASAB1_filtered.pileup": 58,
"CBDAASAB2_filtered.pileup": 42,
"ALQAAAD2_filtered.pileup": 132,
"GARAAAD5_filtered.pileup": 96
}


fmax_coverage = {
"MONAASAB3_filtered.pileup" : 74,
"BAGAASAB1_filtered.pileup" : 72,
"CBDAASAB2_filtered.pileup" : 54,
"ALQAAAD2_filtered.pileup" : 78,
"GARAAAD5_filtered.pileup" :60
}

path_to_gtf  ="./gtfs"

gts_to_parse = [x for x in os.listdir(path_to_gtf) if x.endswith(".gtf")]

done_filtered = [x for x in os.listdir("./piatpositions/resultsfiltered/pi") if x.endswith(".pi")] +\
                [x for x in os.listdir("./piatpositions/resultsfiltered/theta") if x.endswith(".theta")] +\
                [x for x in os.listdir("./piatpositions/resultsfiltered/D") if x.endswith(".D")]


done_filtered = set(done_filtered)

for gtf_file in gts_to_parse:
    gtf_path = f"{path_to_gtf}/{gtf_file}"


    for fpileup_file in fmax_coverage.keys():
        fpileup_path = f"./pianalysis/pileups/{fpileup_file}"

        fcoverage = fmax_coverage[fpileup_file]
        fpool_size= fpool_sizes[fpileup_file]

        for fmeasure in ["pi","theta","D"]:

            print(f"starting {gtf_file} {fpileup_file} {fmeasure}")
            foutput_file = f"{fpileup_file.replace('.pileup','')}_{gtf_file}_measure_{fmeasure}.{fmeasure}".replace("_subset_individual_CpG","")

            if foutput_file in done_filtered:
                print("skipping {}".format(foutput_file))
                continue


            foutput_path = f./piatpositions/resultsfiltered/{fmeasure}/{foutput_file}"

 
            subprocess.run(["perl","./popoolation_1.2.2/Variance-at-position.pl",\
                        "--pool-size",str(fpool_size), "--min-qual", "20", "--min-count","2",\
                        "--min-coverage","5","--max-coverage",str(fcoverage),\
                        "--fastq-type", "sanger", "--pileup",fpileup_path, "--gtf",gtf_path, "--measure",fmeasure,\
                        "--output",foutput_path ])


