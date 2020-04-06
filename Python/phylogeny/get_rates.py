from __future__ import division
import os
import phyphy
import pandas as pd
from Bio import SeqIO

mydir = os.path.expanduser("~/GitHub/Evol_16S/")

def get_rates():
    #hyphy_path = '/usr/local/bin/HYPHYMPI'
    #hyphy_path = "/usr/local/lib/hyphy/"
    #h = HyPhy(build_path = hyphy_path, suppress_log = True, quiet=True)
    h = phyphy.HyPhy(executable   = "HYPHYMPI",
                           mpi_launcher = "mpirun",
                           mpi_options  = "-np 32")

    align = mydir + 'data/phylogeny/HYPHYMP/nmicrobiol201648-s7_clean_noEuks.txt'
    tree = mydir + 'data/phylogeny/HYPHYMP/RAxML_bestTree.T20'
    out = mydir + 'data/phylogeny/HYPHYMP/rates.json'
    r = phyphy.LEISR(hyphy=h, alignment = align, tree = tree, type = "Nucleotide", output = out, model = "GTR", rate_variation = 'Gamma')
    r.run_analysis()
    #p = HyPhyParser(json_path = "sim.json")
    #p.extract_csv(hyphy_homo_csv)

def clean_out():
    IN_path = mydir + 'data/phylogeny/HYPHYMP/out.txt'
    OUT_path = mydir + 'data/phylogeny/HYPHYMP/out_clean.txt'
    OUT = open(OUT_path, 'w')
    header = 'Site\tRate\tCI_025\tCI_975\n'
    OUT.write(header)
    problem_sites = ['75', '85', '86', '100', '101', '103', '140', '223', \
                    '224', '483', '638', '639', '640', '641', '642', '899', \
                    '902', '1071', '1072', '1073', '1074', '1117', '1114', \
                    '1119', '1134', '1144', '1230', '1567', '1574', '1572', \
                    '1573', '1575', '1571', '1581', '1586', '1568', '1582', \
                    '1583', '1584', '1585', '1566', '1587', '1588', '1592', \
                    '1590', '1576', '1570', '1569', '1589', '1591']
    for line in open(IN_path, 'r'):
        line_split = line.strip().split()
        if len(line_split) < 1:
            continue

        if (line_split[0] == '|') and (line_split[1] != 'Site'):
            site = line_split[1]
            rate = line_split[3]
            rate_025 = line_split[5]
            if site in problem_sites:
                rate_975 = line_split[6][1:]
            else:
                rate_975 = line_split[7]
            line = site + '\t' + rate + '\t' + rate_025 + '\t' + rate_975 +  '\n'
            OUT.write(line)
    OUT.close()

def get_nuc_freqs(seq_length = 1670):
    nucs = ['A', 'C', 'G', 'T']
    rates_path = mydir + 'data/phylogeny/HYPHYMP/out_clean.txt'
    rates = pd.read_csv(rates_path, sep = '\t')
    sites = rates.Site.values
    nuc_dict = {}
    for site in sites:
        nuc_dict[site] = {}
        nuc_dict[site]['A'] = 1
        nuc_dict[site]['C'] = 1
        nuc_dict[site]['G'] = 1
        nuc_dict[site]['T'] = 1
    align_path = mydir + 'data/phylogeny/HYPHYMP/nmicrobiol201648-s7_clean_noEuks.txt'
    for record in SeqIO.parse(align_path, "fasta"):
        # start at 1
        for i, site in enumerate(record.seq, 1):
            if (i in sites) and (site != '-') and (site in nucs):
                nuc_dict[i][site] += 1
    nuc_df = pd.DataFrame(nuc_dict).T
    nuc_df.loc[:,:] = nuc_df.loc[:,:].div(nuc_df.sum(axis=1), axis=0)
    out_path = mydir + 'data/phylogeny/HYPHYMP/site_specific_nucs.txt'
    nuc_df['Site'] = nuc_df.index
    nuc_df.to_csv(out_path, sep = '\t', index = False)



def merge_dfs():
    #merge site specific rates and nuc freqs , make sure they match up
    freqs_path = mydir + 'data/phylogeny/HYPHYMP/site_specific_nucs.txt'
    freqs = pd.read_csv(freqs_path, sep = '\t')
    rates_path = mydir + 'data/phylogeny/HYPHYMP/out_clean.txt'
    rates = pd.read_csv(rates_path, sep = '\t')
    merged = rates.merge(freqs, on =  'Site', how='outer')
    merged_clean = merged.loc[~(merged['CI_025'] == 0.0) & ~(merged['CI_975'] >= 4600)]
    out_path = mydir + 'data/phylogeny/HYPHYMP/rates_freqs.txt'
    merged_clean.to_csv(out_path, sep = '\t', index = False)



#get_rates()
#clean_out()
#get_nuc_freqs()
#merge_dfs()
