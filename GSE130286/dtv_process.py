import GEOparse
import pandas as pd
import sqlite3
import csv
import os
import sys
import math
import re
import tempfile
import statistics
import argparse
import pymongo
import gzip
import shutil
import numpy as np
import io
import scanpy as sc
import glob
from anndata import AnnData

accessionID = 'GSE129800'        # GSE ID
experimentType = 'm'     # m or rna      (microarray or RNA-Seq)
norm = 'None'               # tpm/cpm/None  Need to normalize data (usually TPM)
takeLog = 'True'            # True/False    need to take log of data (data is not log normalized yet)
raw = 'False'                # True/False    RNA-Seq raw data?

if accessionID == '' or experimentType == '' or takeLog == '' or norm == '' or raw == '':
    print("System usage: [GEO Accession ID] [m (microarray) or rna (RNA-Seq)] [takelog: True/False] [Normalization: cpm or tpm or none] [raw seq: True/False]")
    quit()

# GEOParse

path='./'+str(accessionID)+'_family.soft.gz'
path_dir='./'

if os.path.exists(path):
    # Load from existing file
    print("-Loading from", path)
    gse = GEOparse.get_GEO(filepath = path, silent=False)

else:
    # Download GSE from GEO and load it
    print('-Downloading', str(accessionID))

    gse = GEOparse.get_GEO(geo=str(accessionID), destdir=path_dir)
    
print("Done processing")


with gzip.open(path_dir+str(accessionID)+'_family.soft.gz', 'rb') as f_in:
    with open(path_dir+str(accessionID)+'_family.soft.txt', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)


def load(filename, con, table):
    with open(filename) as tsvfile:
        reader = csv.reader(tsvfile, dialect='excel-tab')
        header = True
        for row in reader:
            r = row[:]
            if header:
                query = 'CREATE TABLE IF NOT EXISTS %s (%s)' % (table, ', '.join(["`%s` text" % s for s in r[:]]))
            else:
                query = 'INSERT INTO %s VALUES (%s)' % (table, ', '.join(["\"%s\"" % v for v in r[:]]))
            header = False
            con.execute(query)

"""

Takes log 2 of expr

"""

def logtake(val):
    v = float(val)
    if v > 1:
        result = np.log2(v)
    elif v > -1:
        result = (v) - 1
    else:
        result = -2 - np.log2(-1 * v)
    return result

"""

Generates expr.txt file for microarray experiments

"""

def m_analyze(gpl):
    print('Starting m_analyze')
    # Create an empty dataframe with rows and cols
    rownames = []
    colnames = []
    for name, gsm in gse.gsms.items():
        if gsm.metadata['platform_id'][0] != gpl.name:
            continue
        colnames.append(name)
        if len(rownames) == 0:
            rownames = gsm.table['ID_REF'].values.tolist()

    gse_df = pd.DataFrame(index=rownames, columns = colnames)
    for samplekey, sample in gse.gsms.items():
        if sample.metadata['platform_id'][0] != gpl.name:
            continue
        print(samplekey)
        column = sample.table['VALUE'].tolist()
        
        # If values are in log 10 form, need to change to log 2
        #column = (sample.table['VALUE']/math.log(10, 2)).tolist()

        if takeLog == 'True':
            print("Taking log")
            column = [np.math.log(y, 2) if y>0 else -1 for y in column]
            gse_df[samplekey] = column

        if(len(gse_df[samplekey]) != len(column)):
            print('not equal in %s'%samplekey)
            print(len(gse_df[samplekey]))
            print(len(column), type(column))
            filler = -1
            gse_df.loc[:, samplekey] = column + [filler] * (len(gse_df.index) - len(column))
        else: 
            gse_df[samplekey] = column

    gse_df.insert(0, 'ID', rownames)

    gse_df.to_csv(path_dir+str(accessionID)+'gse_withid-%s.txt'%gpl.name, header=True, index=False, sep='\t')

    if os.path.exists(path_dir+str(accessionID)+'-%s.txt'%gpl.name):
        os.remove(path_dir+str(accessionID)+'-%s.txt'%gpl.name)
    print(gpl.table.columns)
    col_symbol = []
    col_title = []
    
    if 'gene_assignment'in gpl.table.columns:
        lst = gpl.table['gene_assignment']
        for l in lst:

            if l != '---':
                l = str(l)
                x = [y.split(' // ') for y in l.split(' /// ')]
                col_title.append(" /// ".join(i[2] if len(i) > 2 else "ERR" for i in x))
                col_symbol.append(" /// ".join(i[1] if len(i) > 2 else "ERR" for i in x))
            else:
                col_symbol.append("---")
                col_title.append("---")
                
    # check it, Definition and symbol can be saved with different name other than the following
    if 'GeneSymbol'in gpl.table.columns:
        gpl.table=gpl.table.rename(columns = {'GeneSymbol':'Symbol'})
    if 'Gene Symbol'in gpl.table.columns:
        gpl.table=gpl.table.rename(columns = {'Gene Symbol':'Symbol'})
    if 'gene symbol'in gpl.table.columns:
        gpl.table=gpl.table.rename(columns = {'gene symbol':'Symbol'})

    if 'geneName'in gpl.table.columns:
        gpl.table=gpl.table.rename(columns = {'geneName':'Definition'})
    if 'Gene Name'in gpl.table.columns:
        gpl.table=gpl.table.rename(columns = {'Gene Name':'Definition'})
    if 'Gene Title'in gpl.table.columns:
        gpl.table=gpl.table.rename(columns = {'Gene Title':'Definition'})
    if 'Gene Description' in gpl.table.columns:
        gpl.table = gpl.table.rename(columns={'Gene Description': 'Definition'})
    if 'Description' in gpl.table.columns:
        gpl.table = gpl.table.rename(columns={'Description': 'Definition'})
    if 'DESCRIPTION' in gpl.table.columns:
        gpl.table = gpl.table.rename(columns={'DESCRIPTION': 'Definition'})
    if 'GENE_SYMBOL' in gpl.table.columns:
        gpl.table = gpl.table.rename(columns={'GENE_SYMBOL': 'Symbol'})
    if 'GENE_SYM' in gpl.table.columns:
        gpl.table = gpl.table.rename(columns={'GENE_SYM': 'Symbol'})
    if 'description' in gpl.table.columns:
        gpl.table = gpl.table.rename(columns={'description': 'Definition'})

    if 'gene_assignment'in gpl.table.columns:
        gpl.table['Symbol'] = col_symbol
        gpl.table['Definition'] = col_title

    
    gpl.table['Symbol'] = gpl.table['ID']
    gpl.table['Definition'] = gpl.table['ID']
    gpl.table.to_csv(path_dir+str(accessionID)+'-%s.txt'%gpl.name, header=True, index=False, sep='\t')

    if os.path.exists('example.db'):
        os.remove('example.db')
    conn = sqlite3.connect('example.db')
    c = conn.cursor()

    load(path_dir+str(accessionID)+'-%s.txt'%gpl.name, c, gpl.name)
    load(path_dir+str(accessionID)+'gse_withid-%s.txt'%gpl.name, c, 'gsed_%s'%gpl.name)

    c.execute("drop TABLE IF EXISTS new_gse")
    c.execute("CREATE TABLE new_gse AS SELECT '' as ProbeID, '' as Name, * FROM gsed_%s where 0"%gpl.name)

    c.execute("""INSERT INTO new_gse SELECT gsed_{0}.ID as ProbeID,
                        {0}.Symbol || ':' || {0}.Definition as Name,
                        gsed_{0}.*
                 FROM {0},gsed_{0} WHERE {0}.ID=gsed_{0}.ID""".format(gpl.name))

    conn.commit()
    conn.close()

    # Generate expr.txt
    con = sqlite3.connect('example.db')
    data_expr=pd.read_sql_query("SELECT * FROM new_gse", con)
    del data_expr['ID']

    data_expr.to_csv(path_dir+str(accessionID)+'-%s-expr.txt'%(gpl.name), header=True, index=False, sep='\t')
    con.close()
    print("Done with m_analyze")

""" file2dict()


"""

def file2dict(filename, id_col, prefix=''):
    table =collections.defaultdict(dict)
    reader = csv.reader(filename, delimiter = '\t', dialect='excel')
    
    columns = next(reader) # Get top row
    id_idx = columns.index(id_col) # Find col of data needed
    ids = []

    for row in reader:
        cid = row[id_idx]
        ids.append(cid)
        for col in columns:
            if col == id_col:
                continue
            table[cid][prefix+col] = row[columns.index(col)]

    return table, ids, [prefix+c for c in columns]

"""

Generates expr file for RNA-seq experiments

"""

def rna_analyze(gpl):
  #  os.rename('./to_process/'+directory+'/'+'raw/', path_dir+'raw/')
    if raw == 'True':
        files = glob.glob(path_dir+'raw/*.txt.gz')
        print('j', files)
        df = None
        for f in files:
            arr = re.sub(".*(GSM[0-9]+).*", "\\1", f)
            df1 = pd.read_csv(f, sep="\t", header=None, index_col=0, names=['Gene', arr])
            if 'no_feature' in df1.index.values:
                df1 = df1.drop('no_feature')
            if 'ambiguous' in df1.index.values:
                df1 = df1.drop('ambiguous')
            if 'too_low_aQual' in df1.index.values:
                df1 = df1.drop('too_low_aQual')
            if 'not_aligned' in df1.index.values:
                df1 = df1.drop('not_aligned')
            if 'alignment_not_unique' in df1.index.values:
                df1 = df1.drop('alignment_not_unique')
            if df is None:
                df = df1
            else:
                if arr in df.columns:
                    df[arr] += df1[arr]
                else:
                    df = df.merge(df1, 'outer', on='Gene')
        df = df.replace(np.NaN, 0)
    
        df.to_csv(path_dir+str(accessionID)+"-counts.txt", sep='\t')

    else:
        df = pd.read_csv(path_dir+str(accessionID)+"-counts.txt", sep='\t')
        

    expr = df.copy(deep=True)
    expr = expr.drop(['ProbeID', 'Name'], axis=1)

    print("Normalizing")
    adata = AnnData(expr.T)
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata, base=2)

    norm_df = pd.DataFrame(adata.X)
    norm_df = norm_df.T
    
    if raw == 'True':
        norm_df.insert(0, 'ProbeID', df.index.values.tolist())
        norm_df.insert(1, 'Name', df.index.values.tolist())
    
    else:
        norm_df.insert(0, 'ProbeID', list(df['ProbeID']))
        norm_df.insert(1, 'Name', list(df['Name']))
        norm_df.columns = list(df.columns)        

    norm_df.to_csv(path_dir+str(accessionID)+'-%s-expr.txt'%(gpl.name), header=True, index=False,sep='\t')
    
    print("Done with rna_analyze")


""" make_idx: Generates idx.txt

"""

def make_idx(gpl):
    print('Starting make_idx')
    expr = path_dir+str(accessionID)+'-%s-expr.txt'%(gpl.name)

    ptr = []
    ids = []
    name = []
    desc = []
    pos = 0

    with open(expr, 'rb') as f:
        for line in f:
            if pos == 0:
                pos += len(line)
            else:
                ptr.append(pos)
                pos += len(line)
                split = line.decode("utf-8").split('\t')
                ids.append(split[0])
                name.append(split[1].split(':')[0])
                desc.append(':'.join(split[1].split(':')[1:]))
        f.close()

    with open(path_dir+str(accessionID)+'-%s-idx.txt'%gpl.name, 'w') as f:
        f.write('ProbeID\tPtr\tName\tDescription\n')
        for i in range(len(ids)):
            f.write('{}\t{}\t{}\t{}\n'.format(ids[i], ptr[i], name[i], desc[i]))
        f.close()
    print("Done with make_idx")

""" survival_ih: generate survival and ih files

"""

def survival_ih(gpl):
    print('Starting survival_ih')
    survival = open(path_dir+str(accessionID)+'-%s-survival.txt'%gpl.name, 'w')
    flname = path_dir+str(accessionID)+'_family.soft.txt'

    GSM_order = []
    sampleID = {} # GSM numbers
    fileptr = {}  # probeid: fileptr

    # array for headers (ih file)
    array_header = []

    # array for clinical header (ih file)
    clinical_header = []

    # list of headers of channel info/columns
    survival_header_list = ["ArrayId", 'time', 'status', 'title']

    # prefix/type of input of the column
    survival_header_pre = ['', '', '', 'c ']

    # declare channel variables
    channel_num = ''
    channel = ''

    # whether keys have been added to the dictionary
    header_added = False

    # initialize survial info dictionary
    survival_info = {}
    survival_info['ArrayId'] = []
    survival_info['time'] = []
    survival_info['status'] = []
    survival_info['title'] = []

    with io.open(flname, 'r', encoding='utf8') as fe:
        for line in fe:
            if '^SAMPLE' in line:
                flag_gpl = True
                id_string = line.strip('^SAMPLE = ').strip('\n')
                platform = gse.gsms[id_string].metadata['platform_id']
                if platform[0] != str(gpl.name):
                    flag_gpl = False
                    continue
                else: True
                if flag_gpl:
                    GSM_order.append(id_string)
                    array_header.append(id_string)
                    sampleID.update({id_string: []})
                    fileptr[id_string] = []
                    survival_info['ArrayId'].append(id_string)

            if "!Sample_series_id" in line and (flag_gpl):
                series_id = line.replace("!Sample_series_id = ", '').strip('\n')
            elif '!Sample_title' in line and (flag_gpl):
                c_header = line.split(' = ')[1].strip('\n')
                clinical_header.append(c_header)
                survival_info['title'].append(c_header)

            if "!Sample_channel_count" in line and (flag_gpl):
                channel_count = line.strip("!Sample_channel_count = ").strip('\n')
                
                while True:
                    ln = next(fe)

                    if "!Sample_treatment_protocol_ch1" in line:
                        print("found")

                    # Break when line does not  contain _ch#
                    if '_ch' not in ln:
                        header_added = True
                        break
                    if 'protocol' in ln:
                        continue
                    if 'label' in ln:
                        continue

                    ch = int(re.search(r'_ch(\d)',ln).group(1))
                    channel = '_ch%d' %ch
                    
                    # Gets the header/key
                    lst = [x.strip() for x in ln.split(' = ')]
                    titles = lst[0].replace("!Sample_",  '')
                    title = '%s (ch%d)' % (titles.split(channel)[0], ch)

                    # If ':'  exists in the second part, the substring before
                    # ':' is header and after is the value
                    if ':' in  lst[1]:
                        titles = [x.strip() for x in lst[1].split(':')]
                        title =  '%s (ch%d)' % (titles[0], ch)
                        lst[1] = titles[1]

                    # If header has not been added, add to dictionary, check
                    # its type (char or num) and append to prefix list
                    if not header_added:
                        survival_header_pre.append('c ')
                        survival_header_list.append(title)
                    if title not in survival_info:
                        survival_info[title] = [lst[1].strip('\n')]
                    else:
                        survival_info[title].append(lst[1].strip('\n'))

    # IH Writing
    ih_file = open(path_dir + str(accessionID) + '-%s-ih.txt' % (gpl.name), 'w')
    ih_header = 'ArrayID' + '\t' + 'ArrayHeader' + '\t' + 'ClinicalHeader'
    ih_file.write(ih_header)
    ih_file.write('\n')

    # Write in ih file
    counter = 0
    print('array_header', array_header)
    print('clinical_header', clinical_header)
    for element in GSM_order:
        ih_file.write('\t'.join([element, array_header[counter].replace(u"\u2019", "'"), clinical_header[counter].replace(u"\u2019", "'")]))
        ih_file.write('\n')
        counter  += 1

    ih_file.close()

    i = 0
    length = len(survival_header_list)
    survival_header = ''
    # Write header
    for i in range(length):
        # No tab after the last header
        if i == length - 1:
            survival_header =  survival_header + survival_header_pre[i] + survival_header_list[i]
        else:
            survival_header = survival_header + survival_header_pre[i] + survival_header_list[i] + '\t'
    survival_header = survival_header + '\n'
    survival.write(survival_header)

    # Write body
    i = 0
    length = len(survival_info['ArrayId'])
    for i in range(length):
        # a line
        aTuple = []
        for t in survival_header_list:
            # If no value is added to the column, append empty string
            if len(survival_info[t]) == 0:
                aTuple.append('')
            else:
                try:
                    aTuple.append(survival_info[t][i])
                except IndexError:
                    aTuple.append('')
                continue

        # write a line
        survival.write('\t'.join([str(s) for s in aTuple]))
        survival.write('\n')
    survival.close()
    print("Done with survival_ih")
                     

""" main(): runs the program

Args:
    accessionID: GEO accession ID
    takeLog: True if need to take log of data

"""

for _, gpl in gse.gpls.items():
    
    if experimentType == 'm':
        m_analyze(gpl)
    else:
        rna_analyze(gpl)
    survival_ih(gpl)
    make_idx(gpl)
