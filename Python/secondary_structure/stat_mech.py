from __future__ import division
import os, re, math, random
import pandas as pd
import subprocess as sp
from random import randint


mydir = os.path.expanduser('~/github/Evol_16S/')


u_m = 2.3
u_mm = 2.3

energy_dict = { ('A', 'U'): -2*u_m,  ('G', 'C'): -3*u_m, \
        ('U', 'A'): -2*u_m, ('C', 'G'): -3*u_m}

bp_keep = [['A','U'], ['U','A'], ['C','G'], ['G','C']]

bad_rnaeval_bp = [['G','A'], ['A','G'],['A','C'], ['C','A'], \
                ['U','U'], ['A','A'], ['C','C'], ['G','G'], \
                ['U','C'], ['C','U']]


def get_size(file_path):
    df = pd.read_csv(file_path, sep = "\t")
    return df.shape[0]


def get_mfe(df):
    """Runs RNAeval for RNA sequence and its structure.
    Returns the mfe value.
    """
    #df_path = mydir + 'data/secondary_structure/no_pseudoknot/SB.bacilli.merged/AB021192.txt
    seq = "".join(df.Base.values)
    struct = "".join(df.Bracket.values)
    rnaeval_input = seq + '\n' + struct
    bad_bp_count = 0
    good_bp_count = 0
    loop_bp_count = 0
    for index, row in df.iterrows():
        if row['Paired_base'] != 0:
            base_1 = row['Base']
            base_2 = df.ix[row['Paired_base'] -1]['Base']
            #print(base_1, base_2)
            if [base_1, base_2] in bp_keep:
                good_bp_count += 1
            else:
                bad_bp_count += 1
        else:
            loop_bp_count += 1
    rnaeval_proc = sp.Popen("RNAeval", shell=True, stdin=sp.PIPE,
                            stdout=sp.PIPE)
    #stdout=FNULL, stderr=subprocess.STDOUT
    #sp.PIPE
    rnaeval_stdout = rnaeval_proc.communicate(rnaeval_input.encode())
    # RNAeval returns the following string:
    # 'sequence\nstructure (mfe_value)\n'
    mfe = float(rnaeval_stdout[0].decode().strip().split()[2][1:][:-1])
    percent_bad_bp = bad_bp_count / (bad_bp_count + good_bp_count)
    percent_stem_all = (bad_bp_count + good_bp_count) / df.shape[0]
    percent_stem_good = good_bp_count / df.shape[0]
    return (mfe, percent_bad_bp, percent_stem_all, percent_stem_good)


def get_mfe_all():
    directory = mydir + 'data/secondary_structure/no_pseudoknot/bacilli/SB.bacilli.merged/'
    OUT_path = mydir + 'data/secondary_structure/no_pseudoknot/mfe.txt'
    OUT = open(OUT_path, 'w')
    header = 'Sample\tSize\tMFE\tPercent_bad\tPercent_stem_all\tPercent_stem_good\n'
    OUT.write(header)
    for x in os.listdir(directory):
        if x.endswith(".txt"):
            sample_name = x.rsplit('.', 1)[0]
            df = pd.read_csv(os.path.join(directory, x), sep = "\t")
            mfe = get_mfe(df)
            size = get_size(os.path.join(directory, x))
            line = sample_name + '\t' + str(size) + '\t' + str(mfe[0]) + '\t' + \
                    str(mfe[1]) + '\t' + str(mfe[2]) + '\t' + str(mfe[3]) +  '\n'
            OUT.write(line)
    OUT.close()


def get_local_rna(df, site):
    plus_one_type = df.ix[site+1,'Bracket']
    # if plus_one_type == '(' matching bp is downstream
    # if plus_one_type == ')' matching bp is upstream
    #plus_one_paired_index = (df.ix[site+1, 'Paired_base']) - 1
    print(df.ix[289,])
    print(site+1)
    print(plus_one_type)
    print(plus_one_paired_index)


def get_mutant(df):
    # this algorithm is technically JC69 corrected
    nucs = ['A', 'C', 'G', 'T']
    #site = randint(0, df.shape[0])
    site = 304
    old_base = df.ix[site,'Base']
    if old_base in nucs:
        nucs.remove(old_base)
    new_base = random.choice(nucs)
    df.ix[site,'Base'] = new_base
    if site > 0 and site < df.shape[0]:
        if df.ix[site,'Bracket'] != '.':
            paired_base = (df.ix[site, 'Paired_base'])
            df.ix[site,'Bracket'] = '.'
            df.ix[paired_base-1,'Bracket'] = '.'
            df.ix[site,'Paired_base'] = 0
            df.ix[paired_base-1,'Paired_base'] = 0
            #df.ix[paired_base - 1,'Paired_base'] = df.ix[site,'Index']
        else:
            if df.ix[site+1,'Bracket'] == ')' or df.ix[site+1,'Bracket'] == '(':
                #get_local_rna(df, site)
                plus_one_type = df.ix[site+1,'Bracket']
                # if plus_one_type == '(' matching bp is downstream
                # if plus_one_type == ')' matching bp is upstream
                plus_one_paired_index = (df.ix[site+1, 'Paired_base']) - 1
                site_paired_bracket_index = plus_one_paired_index + 1
                site_paired_bracket = df.ix[site_paired_bracket_index,'Bracket']
                #if site_paired_bracket
                #print(site_paired_bracket)
                #print(df.ix[plus_one_paired_index:site+2,])


                #df.ix[site,'Bracket'] = df.ix[site+1,'Bracket']
                #print(plus_one_type)
                #print(paired_base_index)
                #print(df.ix[site:site+20,])



                #print(df.ix[site,])
                #df.ix[site,'Paired_base'] = df.ix[paired_base_index_plus,'Index']
                #print(df.ix[site,])
                #print(df.ix[paired_base_index_plus ,])
                #print(df.ix[paired_base_index_plus + 1,])
                #df.ix[paired_base_index_plus ,'Paired_base'] = df.ix[site,'Index']
                #df.ix[paired_base_index_plus ,'Bracket'] = df.ix[paired_base_index_plus-1,'Bracket']
                #print(site)
                #print('plus_bracket')
                #print(df.ix[site-20:site+20,])
                #print(df.ix[site,])
                #print(df.ix[paired_base_index_plus,])

            elif df.ix[site-1,'Bracket'] == ')' or df.ix[site-1,'Bracket'] == '(':
                #print(df.ix[site-20:site+20,])
                print(site)
                print('minus_bracket')
                # -1 to account for zero-based indexing
                paired_base_index_minus = (df.ix[site-1, 'Paired_base']) - 1
                df.ix[site,'Bracket'] = df.ix[site-1,'Bracket']
                df.ix[site,'Paired_base'] = df.ix[paired_base_index_minus-1,'Index']
                df.ix[paired_base_index_minus - 1,'Paired_base'] = df.ix[site,'Index']
                df.ix[paired_base_index_minus - 1,'Bracket'] = df.ix[paired_base_index_minus,'Bracket']
                #print(df.ix[site,])
                #print(df.ix[paired_base_index_minus,])

    return df


def get_mutant_rnafold(df):
    nucs = ['A', 'C', 'G', 'T']
    #site = randint(0, df.shape[0])
    site = 304
    old_base = df.ix[site,'Base']
    if old_base in nucs:
        nucs.remove(old_base)
    new_base = random.choice(nucs)
    df.ix[site,'Base'] = new_base
    if df.ix[site,'Bracket'] != '.':
        paired_base = (df.ix[site, 'Paired_base'])
        df.ix[site,'Bracket'] = '.'
        df.ix[paired_base-1,'Bracket'] = '.'
        df.ix[site,'Paired_base'] = 0
        df.ix[paired_base-1,'Paired_base'] = 0
    else:
        bracket = df.Bracket.values
        paired_base = df.Paired_base.values
        if df.ix[site+1,'Bracket'] == ')' or df.ix[site+1,'Bracket'] == '(':
            plus_one_type = df.ix[site+1,'Bracket']
            paired_base_x = paired_base[site + 1]
            sites_subset = [site+1]
            for x in range(site+2,len(bracket)):
                if bracket[x] == plus_one_type and ( paired_base[x] == paired_base_x +1  or paired_base[x] == paired_base_x -1 ):
                    sites_subset.append(x)
                    paired_base_x = paired_base[x]
                else:
                    break
            df_sub = df.ix[sites_subset,]
            start = df_sub.ix[df_sub['Paired_base'].argmin(),'Paired_base'] - 1
            stop = df_sub.ix[df_sub['Index'].argmax(),'Index'] - 1
            df_sub_all = df.iloc[start:stop,]
            seq = "".join(df_sub_all.Base.tolist()) + '\n'

            rnaf = sp.Popen(["RNAfold","--noPS"],
                            stdin=sp.PIPE,
                            stdout=sp.PIPE,
                            stderr=sp.PIPE,
                            # Universal Newlines effectively allows string IO.
                            universal_newlines=True)
            foldout = rnaf.communicate(seq)
            #print(df_sub_all)
            #print(foldout)
            #foldout = rnaf.communicate(seq)
            #print(foldout)

            #exp = re.compile(r'\(\s*(-?\d+\.\d+)\)')

            #s = sp.check_output(['RNAfold', '--noPS',],
			#												input="{}\n@\n".format(seq),
			#												universal_newlines=True)

            #m = exp.search(s.split('\n')[1])
            #print(s)
            #f = float(m.groups(1)[0])
            #print(f)
            #print(df_sub_all)
            #get_mfe(df_sub_all)
            # next, make sure there are no extra brackets
            # if there are, add those in, repeat


def test_mutant():
    df_path =  mydir + 'data/secondary_structure/no_pseudoknot/bacilli/SB.bacilli.merged/Ab.anurly3.txt'
    OUT_path = mydir + 'data/secondary_structure/no_pseudoknot/delta_delta_G.txt'
    OUT = open(OUT_path, 'w')
    header = 'Count\tDelta_delta_G\n'
    OUT.write(header)

    df = pd.read_csv(df_path, sep = "\t", index_col = 0)
    mfe = get_mfe(df)
    delta_G_anc = mfe[0]
    count = 0
    while count < 1000:
        df_mut = get_mutant(df)
        mfe_mut = get_mfe(df_mut)
        delta_G_mut = mfe_mut[0]
        if delta_G_mut < 0:
            print(count)
            count += 1
            line = str(count) + '\t' + str(delta_G_mut-delta_G_anc) + '\n'
            OUT.write(line)
    OUT.close()

#test_mutant()
df_path =  mydir + 'data/secondary_structure/no_pseudoknot/bacilli/SB.bacilli.merged/Ab.anurly3.txt'
df = pd.read_csv(df_path, sep = "\t", index_col = 0)
#print(df.ix[304-5:304+7,])
#print(df.ix[280:296,])
#print(df.ix[287:310,])
get_mutant(df)

#print(new.ix[304-5:304+7,])
#print(new.ix[280:296,])

#get_mfe(new)
