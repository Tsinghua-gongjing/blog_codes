import pandas as pd 
import numpy as np
import MySQLdb
import sqlalchemy
from nested_dict import nested_dict
import json
import subprocess,sys,os
from multiprocess import Pool
import traceback
import gj
from math import ceil
import math
from collections import OrderedDict
import logging

def fetch_entry_data(mode='rise_id', id='RISE0256681', species='hs', save_dir='/Share/home/zhangqf5/gongjing/DNA-RNA-Protein-interaction-correlation-12-18/results/overlap/duplex_test', show_link_anno='derived', host="bio04",user="zhangqf5",passwd="123456",db="RISE",port=50033, check=0):

    """
	creat conf: circos.conf, zooms.conf, highlight, links for circos plot
	   mode: query id type, also table column queried in table
	   id: query id
	   species: query species
	   save_dir: output dir of circos plot
	   show_link_anno: link type => derived: all links where the two entry genes involved; source: link of the query entry(duplex)
    """

    if check:
    	svg = save_dir+'/'+id+'.convert.svg'
    	resize_html = save_dir+'/'+id+'.convert.resize.html'
    	if os.path.isfile(svg) and os.path.isfile(resize_html):
    		return 

    id_log_fn = save_dir+'/'+id+'.log'
    LOG = open(id_log_fn,'w')
    sys.stdout = LOG
    sys.stderr = LOG

    gj.printFuncRun("[start fetch query data]")
    gj.printFuncArgs()

    data_dict = nested_dict(3, list)

    # connect local database, get the query cursor
    # db = MySQLdb.connect("localhost","root","","RISE")
    db = MySQLdb.connect(host=host, user=user, passwd=passwd, db=db, port=port)
    cur = db.cursor()
    gj.printFuncRun("[success connect the database]")

    table_dict = {'hs':['human_editing','human_modification','human_clipdb','human_pancan','human_ucscSnp142'],
                  'mm':['mouse_editing','mouse_modification','mouse_clipdb','mouse_ucscSnp142'],
                  'yeast':['human_editing']}
    table_type_dict = {'human_DARNED':'modification','human_RADAR':'modification',
                       'human_RMBaseNm':'modification','human_RMBaseOthers':'modification','human_RMBasePseudoU':'modification','human_RMBasem5C':'modification','human_RMBasem6A':'modification',
                       'human_modification':'modification','mouse_modification':'modification',
                       'human_clipdb':'RBP','human_pancan':'SNV','human_ucscSnp142':'SNV',
                       'mouse_DARNED':'modification','mouse_RADAR':'modification',
                       'mouse_PseudoU':'modification','mouse_RMBasem5C':'modification','mouse_RMBasem6A':'modification',
                       'human_editing':'modification','mouse_editing':'modification',
                       'mouse_clipdb':'RBP','mouse_ucscSnp142':'SNV'}
    category_dict = {'PARIS':'GlobalAnalysis','SPLASH':'GlobalAnalysis','MARIO':'GlobalAnalysis','LIGRseq':'GlobalAnalysis',
                     'CLASH':'TargetedAnalysis','RAP':'TargetedAnalysis','RIAseq':'TargetedAnalysis',
                     'NPInterv3.0':'database','RAIDv2.0':'database','RAINv1.0':'database','benchmarkdata':'benchmarkdata'}
    table_query_ls = table_dict[species]
    biomart_gene_table_dict = {'hs':'biomart_hg38', 'mm':'biomart_mm10'}

    ensembl_gene_id_ls = []
    ensembl_gene_name_ls = []
    ensembl_gene_abstract_type_dict = {}
    savefn_zoom = save_dir+'/'+id+'.zooms.conf'

    # get 2 gene_ids for rise_id, dump RRI duplex region into zooms.conf
    if mode == "rise_id":
    	cur.execute("SELECT * FROM %s WHERE rise_id = '%s'"%('rri', id))
    	for row in cur.fetchall():
        	ensembl_gene_id_ls.append(row[11])
        	ensembl_gene_id_ls.append(row[13])
        	ensembl_gene_name_ls.append(row[12])
        	ensembl_gene_name_ls.append(row[14])
        	ensembl_gene_abstract_type_dict[row[11]] = row[17]
        	ensembl_gene_abstract_type_dict[row[13]] = row[18]
        	zooms_print(row, out=savefn_zoom, extend=100, species='hs', scale = 1000)

    if mode == "gene_id":
    	ensembl_gene_id_ls = [id,id]

    #"""
    cur.execute("SELECT * FROM %s WHERE gene_id = '%s' "%('gencode_biomart_hg38_mm10_yeast',id))
    for row in cur.fetchall():
        ensembl_gene_name_ls.append(row[2])
    #"""
    genename_ls_str = '"'+'","'.join(ensembl_gene_name_ls)+'"'
    #"""
    cur.execute("SELECT * FROM %s WHERE gene_name IN (%s) "%('gencode_biomart_hg38_mm10_yeast',genename_ls_str))
    for row in cur.fetchall():
        ensembl_gene_id_ls.append(row[1])
    #"""

    id_ls_str = '"'+'","'.join(ensembl_gene_id_ls)+'"'
    
    print "[interaction gene ids]", ensembl_gene_id_ls, id_ls_str, genename_ls_str

    # get rise entry where the 2/1 gene_ids involves interaction, dump into RRI entry interaction links(links), label duplex genes(not all genes) of query entry
    savefn_link = save_dir+'/'+id+'.link'
    savefn_genes = save_dir+'/'+id+'.genes'
    OUT_LINK = open(savefn_link,'w')
    species_dict = {'hs':'human', 'mm':'mouse'}
    cur.execute("SELECT * FROM %s WHERE (ensembl_gene_id1 IN (%s) OR ensembl_gene_id2 IN (%s)) AND species='%s' "%('rri', id_ls_str, id_ls_str, species_dict[species]))
    gene_link_rise_id_ls = []
    gene_link_gene_id_ls = []
    gene_link_gene_name_ls = []
    gene_link_gene_name_uniq = {}
    for row in cur.fetchall():
    	row = list(row)
        row_uniq1 = ''.join([row[12],row[25],row[14],row[26],row[19],row[20],row[21]])
        row_uniq2 = ''.join([row[14],row[26],row[12],row[25],row[19],row[20],row[21]])
        if gene_link_gene_name_uniq.has_key(row_uniq1) or gene_link_gene_name_uniq.has_key(row_uniq2):
            continue
        else:
            gene_link_gene_name_uniq[row_uniq1] = 1
            gene_link_gene_name_uniq[row_uniq2] = 1
        if row[28] == "." or row[31] == ".":  # filtering entry with no coordinate information, entry of database will not show
        	print "[filtering: no coordinate]",row
        	continue
        if row[28] == "." and row[29] != ".":
        	row[28] = str(int(row[29])-1)
        if row[28] != "." and row[29] == ".":
        	row[29] = str(int(row[28])+1)
        if row[31] == "." and row[32] != ".":
        	row[31] = str(int(row[32])-1)
        if row[31] != "." and row[32] == ".":
        	row[32] = str(int(row[31])+1)
        if mode == "rise_id" and row[7] == id:
        	duplex = "source"
        	with open(savefn_genes,'w') as OUT_GENES:
        		print >>OUT_GENES,' '.join(row[1:4]).replace('chr',species)+' '+row[12]  # interaction genes(not all): chr_stem start_stem end_stem gene_name; not genes' coordinate
        		print >>OUT_GENES,' '.join(row[4:7]).replace('chr',species)+' '+row[14]
        else:
        	duplex = "derived"
        #duplex = "source" if row[7] == id else "derived"
        if species == "yeast":
            print >>OUT_LINK,' '.join('yeast'+i if n==0 or n==3 else i for n,i in enumerate(row[27:33])).replace('chr',species)+' category=%s,method=%s,duplex=%s,gene_id1=%s,gene_name1=%s,gene_type1=%s,gene_id2=%s,gene_name2=%s,gene_type2=%s,id=%s'%(category_dict[row[19]],row[19],duplex,row[11],row[12],row[17],row[13],row[14],row[18],'|'.join(map(str,row[1:])).replace(',','_'))
        else:
            print >>OUT_LINK,' '.join(row[27:33]).replace('chr',species)+' category=%s,method=%s,duplex=%s,gene_id1=%s,gene_name1=%s,gene_type1=%s,gene_id2=%s,gene_name2=%s,gene_type2=%s,id=%s'%(category_dict[row[19]],row[19],duplex,row[11],row[12],row[17],row[13],row[14],row[18],'|'.join(map(str,row[1:])).replace(',','_').replace(' ',''))
        gene_link_rise_id_ls.append(row[7])
        gene_link_gene_id_ls.append(row[11])
        gene_link_gene_id_ls.append(row[13])
        gene_link_gene_name_ls.append(row[12])
        gene_link_gene_name_ls.append(row[14])
        ensembl_gene_abstract_type_dict[row[11]] = row[17]
        ensembl_gene_abstract_type_dict[row[13]] = row[18]
    OUT_LINK.close()
    gene_link_gene_id_ls = list(set(gene_link_gene_id_ls))  # all uniq interacted gene ids
    gene_link_gene_name_ls = list(set(gene_link_gene_name_ls))
    print "[gene_link_gene_id_ls]", gene_link_gene_id_ls
    print "[gene_link_gene_name_ls]", gene_link_gene_name_ls
    print "[rri statistic of %s] involved entry: %s, involved genes: %s"%(id,len(gene_link_rise_id_ls),len(gene_link_gene_id_ls))
    print 

    # get all linked gene's info, dump into .linkgenes
    savefn_linkgenes = save_dir+'/'+id+'.linkgenes'
    linkgenes_info_ls = []
    linkgenes_id_str = '"'+'","'.join(gene_link_gene_id_ls)+'"'
    linkgenes_name_str = '"'+'","'.join(gene_link_gene_name_ls)+'"'

    #cur.execute("SELECT * FROM %s WHERE gene_id IN (%s)"%('biomart_hg38_mm10_longest', linkgenes_id_str))
    #cur.execute("SELECT * FROM %s WHERE gene_id IN (%s)"%('gencode_hg38_mm10', linkgenes_id_str))
    #cur.execute("SELECT * FROM %s WHERE gene_id IN (%s)"%('gencode_biomart_hg38_mm10_yeast', linkgenes_id_str))
    cur.execute("SELECT * FROM %s WHERE gene_id IN (%s)"%('gencode_biomart_hg38_mm10_yeast', linkgenes_id_str))
    for row in cur.fetchall():
    	linkgenes_info_ls.append(row)
    #print linkgenes_info_ls
    gj.printFuncRun("[fetch all interacted gene's info, prepare processing]")
    if len(linkgenes_info_ls) == 0:
    	gj.printFuncRun("Cannot processing due to len(linkgenes_info_ls)=0")
    	return
    chromosomes_str, chromosomes_scale_str = highlight_print(linkgenes_info_ls, out=savefn_linkgenes, ensembl_gene_abstract_type_dict=ensembl_gene_abstract_type_dict, species=species, mode='linkgenes', duplex_genes=ensembl_gene_id_ls, gene_link_gene_id_ls=gene_link_gene_id_ls)

    # get all linked gene's structure, dump into .genestructure
    savefn_genestructure = save_dir+'/'+id+'.genestructure'
    #genestructure_info_ls = []
    #cur.execute("SELECT * FROM %s WHERE gene_id IN (%s) "%('biomart_hg38_mm10_longest', linkgenes_id_str))
    #for row in cur.fetchall():
    #	genestructure_info_ls.append(row)
    genestructure_str, genestructure_scale_str = genestructure_print(rows=linkgenes_info_ls, out=savefn_genestructure, species=species)

    """
    # get gene_ids info, dump into highlight
    savefn_highlight = save_dir+'/'+id+'.highlight'
    gene_info_ls = []
    cur.execute("SELECT * FROM %s WHERE gene_id IN (%s) "%(biomart_gene_table_dict[species], id_ls_str))
    for row in cur.fetchall():
    	gene_info_ls.append(row)
    chromosomes_str_highlight = highlight_print(gene_info_ls+linkgenes_info_ls, out=savefn_highlight, ensembl_gene_abstract_type_dict=ensembl_gene_abstract_type_dict, species='hs', duplex_genes=ensembl_gene_id_ls)
    """

    gene_link_rise_id_1_str = '"'+'","'.join([i+'_1' for i in gene_link_rise_id_ls])+'"'
    gene_link_rise_id_2_str = '"'+'","'.join([i+'_2' for i in gene_link_rise_id_ls])+'"'
    rise_id_str = '"'+'","'.join([id+'_1',id+'_2'])+'"'
    if show_link_anno == "derived":
    	stem_rise_id_str = gene_link_rise_id_1_str+','+gene_link_rise_id_2_str
    if show_link_anno == "source":
    	stem_rise_id_str = rise_id_str
    print "[fetch_entry_data: RISE stem ids will be query in all tables]",stem_rise_id_str

    # dump query annotations into annotations
    savefn_anno = save_dir+'/'+id+'.anno'
    gj.printFuncRun("[query table ls]"),table_query_ls
    OUT = open(savefn_anno, 'w')
    for table in table_query_ls:
    	cur.execute("SELECT * FROM %s WHERE rise_id IN (%s)"%(table, stem_rise_id_str))
    	print "query table: %s"%(table)
    	n = 0
    	for m,row in enumerate(cur.fetchall()):
    		data_dict[id][2][table].append(row)
    		n += 1
    		#print >>OUT, '\t'.join(map(str,row))
    		row_print(row, OUT, table, table_type_dict, species)
    	print "found entry: %s"%(n)

    db.close()
    OUT.close()
    gj.printFuncRun('[dump annotation into] %s'%(savefn_anno))

    #print json.dumps(data_dict, indent=4)
    print gj.printFuncRun("fetch_entry_data: finished fetch query data")


    save_conf = '/Share/home/zhangqf5/gongjing/software/circos-tutorials-0.67/tutorials/5/test/rise_conf/%s.conf'%(id)
    if len(gene_link_gene_id_ls) > 150:
        create_circos_conf(save_conf=save_conf, species=species, chromosome_str=None, chromosomes_scale_str=None, zooms_conf=None, highlight=savefn_linkgenes, genestructure=savefn_genestructure, anno=savefn_anno, genes=savefn_linkgenes+'.txt', link=savefn_link, scale=0)
    else:
        create_circos_conf(save_conf=save_conf, species=species, chromosome_str=chromosomes_str, chromosomes_scale_str=chromosomes_scale_str, zooms_conf=savefn_genestructure+'.introns', highlight=savefn_linkgenes, genestructure=savefn_genestructure, anno=savefn_anno, genes=savefn_linkgenes+'.txt', link=savefn_link)

    circos_plot(conf=save_conf,savefn_dir=save_dir,savefn=id+'.png')

    #"""
    if len(gene_link_gene_id_ls) <= 150:
        id_all_track_coordinate_convert(id=id, species=species)
        convert_html = save_dir+'/'+id+'.convert.html'
        parse_html(html=convert_html)
    else:
        convert_html = save_dir+'/'+id+'.html'
        parse_html(html=convert_html, savefn_html=save_dir+'/'+id+'.convert.resize.html')
    #"""

    LOG.close()
    sys.stdout = sys.__stdout__

    return save_dir+'/'+id+'.convert.resize.html'

def circos_plot(conf=None,savefn_dir=None,savefn=None, run_dir='/Share/home/zhangqf5/gongjing/software/circos-tutorials-0.67/tutorials/5/test/rise_conf/'):
	gj.printFuncRun('circos_plot')
	gj.printFuncArgs()
	subprocess.call(["cd %s"%(run_dir)],shell=True)
	subprocess.call(["/Share/home/zhangqf5/gongjing/software/circos-0.69-4/bin/circos -conf %s -outputdir %s -outputfile %s -debug_group link"%(conf, savefn_dir, savefn)],shell=True)
	#subprocess.call(["circos -conf %s -outputdir %s -outputfile %s -debug_group link"%(conf, savefn_dir, savefn)],shell=True)
	#subprocess.call(["circos -conf %s -outputdir %s -outputfile %s -debug "%(conf, savefn_dir, savefn)],shell=True)
	subprocess.call(["cd -"],shell=True)
	gj.printFuncRun('circos_plot')
	print 

def row_print(row, out, table, table_type_dict, species='hs'):
    if table in ['human_clipdb','mouse_clipdb']:
        print >>out,' '.join(map(str,row[11:14])).replace('chr',species)+' 1 anno_category=%s,anno_source=%s,anno_type=%s,anno_gene_info=%s,anno_full=%s,id=%s'%(table_type_dict[table], 'CLIPDB',    table_type_dict[table], '|'.join(map(str,row[1:11])).replace(' ',''), '|'.join(map(lambda x:str(x).replace(',','&'),row[11:])).replace(' ',''), '|'.join(map(lambda x:str(x).replace(',','&'),row[1:])).replace(' ','') )
    if table in ['human_editing','mouse_editing']:
        print >>out,' '.join(map(str,row[11:14])).replace('chr',species)+' 1 anno_category=%s,anno_source=%s,anno_type=%s,anno_gene_info=%s,anno_full=%s,id=%s'%(table_type_dict[table], row[18],     'editing',             '|'.join(map(str,row[1:11])).replace(' ',''), '|'.join(map(lambda x:str(x).replace(',','&'),row[11:])).replace(' ',''),  '|'.join(map(lambda x:str(x).replace(',','&'),row[1:])).replace(' ','') )
    if table in ['human_ucscSnp142','mouse_ucscSnp142']:
        print >>out,' '.join(map(str,row[11:14])).replace('chr',species)+' 1 anno_category=%s,anno_source=%s,anno_type=%s,anno_gene_info=%s,anno_full=%s,id=%s'%(table_type_dict[table], 'dbSNP',     'SNP',                  '|'.join(map(str,row[1:11])).replace(' ',''), '|'.join(map(lambda x:str(x).replace(',','&'),row[11:])).replace(' ',''),  '|'.join(map(lambda x:str(x).replace(',','&'),row[1:])).replace(' ','') )
    if table in ['human_modification','mouse_modification']:
        print >>out,' '.join(map(str,row[11:14])).replace('chr',species)+' 1 anno_category=%s,anno_source=%s,anno_type=%s,anno_gene_info=%s,anno_full=%s,id=%s'%(table_type_dict[table], 'RMBase',    row[18].replace(',','_'),                '|'.join(map(str,row[1:11])).replace(' ',''), '|'.join(map(lambda x:str(x).replace(',','&'),row[11:])).replace(' ',''),  '|'.join(map(lambda x:str(x).replace(',','&'),row[1:])).replace(' ','') )
    if table in ['human_pancan']:
        print >>out,' '.join(map(str,row[11:14])).replace('chr',species)+' 1 anno_category=%s,anno_source=%s,anno_type=%s,anno_gene_info=%s,anno_full=%s,id=%s'%(table_type_dict[table], 'PanCancer', 'mutation',             '|'.join(map(str,row[1:11])).replace(' ',''), '|'.join(map(lambda x:str(x).replace(',','&'),row[11:])).replace(' ','').replace('[','(').replace(']',')'), '|'.join(map(lambda x:str(x).replace(',','&'),row[1:])).replace(' ','').replace('[','(').replace(']',')') )  # cannot has symbol: []
	#print >>out,' '.join(map(str,row[11:14])).replace('chr',species)+' 1 anno_type=%s,anno_source=%s,anno_gene_info=%s,anno_full=%s,id=%s'%(table_type_dict[table],table,'|'.join(map(str,row[1:11])),'|'.join(map(lambda x:str(x).replace(',','&'),row[11:])).replace(' ',''),'|'.join(map(lambda x:str(x).replace(',','&'),row[1:10])).replace(' ',''))

def zooms_print(row, out, extend=50, species='hs', scale = 1000000):
	stem_left = [row[1].replace('chr',species),int(row[2])-extend,int(row[3])+extend,scale]
	stem_right = [row[4].replace('chr',species),int(row[5])-extend,int(row[6])+extend,scale]
	zooms_str = "<zooms> \n <zoom> \n chr = %s \n start = %sb \n end = %sb \n scale = %s \n </zoom> \n <zoom> \n chr = %s \n start = %sb \n end = %sb \n scale = %s \n </zoom> \n </zooms>"%tuple(stem_left+stem_right)
	with open(out,'w') as OUT:
		print >>OUT,zooms_str

def genestructure_print(rows, out, species='hs', unit=1000000.0):
	gj.printFuncRun('genestructure_print')
	genestructure_str_ls = []
	genestructure_len_dict = nested_dict(2, list)
	genestructure_info_dict = nested_dict(2, list)
	with open(out, 'w') as OUT:
		for n,row in enumerate(rows):
			print >>OUT,"%s%s %s %s gene_id=%s,gene_name=%s,gene_type=%s,tx_id=%s,sequence_type=%s,id=%s"%tuple([species,row[12]]+list(row[13:15])+[row[1],row[2],row[10],row[6],row[11]]+['|'.join(map(str, row[1:]))])
			genestructure_str_ls.append("%s%s[a%s]:%s-%s"%tuple([species]+[row[12],n,row[13]/unit,row[14]/unit]))
			genestructure_len_dict['a%s'%(n)]['len'] = abs(int(row[13])-int(row[14]))
			genestructure_len_dict['a%s'%(n)]['log10len'] = np.log2(abs(int(row[13])-int(row[14])))
			if abs(int(row[13])-int(row[14])) == 0:
				print "[genestructure_print: sequence start=end]",row
			genestructure_info_dict[row[1]]['start'] = [int(row[7])]
			genestructure_info_dict[row[1]]['end'] = [int(row[8])]
			genestructure_info_dict[row[1]]['regions'].append([int(row[13]),int(row[14])])
			genestructure_info_dict[row[1]]['chr'] = [row[12]]
			genestructure_info_dict[row[1]]['len'] = [abs(int(row[8])-int(row[7]))]
			genestructure_info_dict[row[1]]['gene_name'] = [row[2]]
	genestructure_str = ";".join(genestructure_str_ls)
	print '[genestructure_print: genestructure_str]',genestructure_str
	#print genestructure_len_dict
	#print [j['log10len'] for i,j in genestructure_len_dict.items()]
	#print sum([j['log10len'] for i,j in genestructure_len_dict.items()])
	log10len_all = sum([j['log10len'] for i,j in genestructure_len_dict.items()])
	#print "log10len_all: %s"%(log10len_all)
	for i,j in genestructure_len_dict.items():
		genestructure_len_dict[i]['degree'] = j['log10len']/log10len_all
	genestructure_scale_str = ','.join([i+':'+str(j['degree'])+'rn' if j['degree'] !=0 else i+':'+str(0.000001)+'rn' for i,j in genestructure_len_dict.items()])
	print '[genestructure_print: genestructure_scale_str, by degree]',genestructure_scale_str
	for i,j in genestructure_info_dict.items():
		regions_sorted = sorted(j['regions'])
		genestructure_info_dict[i]['regions'] = regions_sorted
		if len(regions_sorted) == 1 and regions_sorted[0][0] == j['start'][0] and regions_sorted[0][1] == j['end'][0]:
			genestructure_info_dict[i]['introns'] = ['NoIntrons']
			pass
		else:
			regions_intron = []
			for n,k in enumerate(regions_sorted):
				if n == len(regions_sorted)-1:
					break
				else:
					intron_start = k[1]+1
					intron_end = regions_sorted[n+1][0]-1
					if k[1]+1 != regions_sorted[n+1][0]:
						regions_intron.append([intron_start,intron_end])
			if len(regions_intron) > 0:
				genestructure_info_dict[i]['introns'] = regions_intron
			else:
				genestructure_info_dict[i]['introns'] = ['NoIntronsHasUTRCDS']
	#print genestructure_info_dict
	INTRON_TXT = open(out+'.introns.txt','w')
	with open(out+'.introns','w') as INTRON:
		print >>INTRON,"<zooms>"
		for i,j in genestructure_info_dict.items():
			print '[genestructure_print: i-j-1, start check]',i,j
			NoIntronsHasUTRCDS = 0
			if j['introns'][0] == "NoIntrons": # no introns, no scale
				print >>INTRON_TXT,"%s%s %s %s gene_id=%s,type=non_intron,gene_name=%s"%(species,j['chr'][0],j['start'][0],j['end'][0],i,genestructure_info_dict[i]['gene_name'][0])
				continue
			if j['introns'][0] == "NoIntronsHasUTRCDS":
				NoIntronsHasUTRCDS = 1
			print '[genestructure_print: i-j-2, has intons ]',i,j
			introns_len_ls = [p[1]-p[0] for p in j['introns']] if NoIntronsHasUTRCDS != 1 else []
			regions_len_ls = [q[1]-q[0] for q in j['regions']]
			all_len_ls = introns_len_ls + regions_len_ls if NoIntronsHasUTRCDS != 1 else regions_len_ls
			all_len_ratio = [l/float(j['len'][0]) for l in all_len_ls]
			all_log2len_ls = [np.log2(l) for l in all_len_ls]
			all_log2len_sum = sum(all_log2len_ls)
			all_log2len_ratio = [l/float(all_log2len_sum) for l in all_log2len_ls]
			all_scale_ls = [r2/r1 for r1,r2 in zip(all_len_ratio,all_log2len_ratio)]
			for n1,k in enumerate(j['introns']):
				if NoIntronsHasUTRCDS: continue
				if (k[1]-k[0])/float(j['len'][0]) > 0.5:
					scale = 0.02/((k[1]-k[0])/float(j['len'][0]))
				elif (k[1]-k[0])/float(j['len'][0]) > 0.3:
					scale = 0.05/((k[1]-k[0])/float(j['len'][0]))
				else:
					scale = 0.2
				#scale = 0.05/((k[1]-k[0])/float(j['len'][0])) if (k[1]-k[0])/float(j['len'][0]) > 0.3 else 0.2
				print >>INTRON,"<zoom>\n chr = %s%s \n start = %sb \n end = %sb \n scale = %s \n </zoom>"%(species, j['chr'][0],k[0],k[1],all_scale_ls[n1])
				print >>INTRON_TXT,"%s%s %s %s gene_id=%s,type=intron,gene_name=%s"%(species,j['chr'][0],k[0],k[1],i,genestructure_info_dict[i]['gene_name'][0])
			if NoIntronsHasUTRCDS: n1 = 0
			for n2,m in enumerate(j['regions']):
				print >>INTRON,"<zoom>\n chr = %s%s \n start = %sb \n end = %sb \n scale = %s \n </zoom>"%(species, j['chr'][0],m[0],m[1],all_scale_ls[n1+n2])
				print >>INTRON_TXT,"%s%s %s %s gene_id=%s,type=non_intron,gene_name=%s"%(species,j['chr'][0],m[0],m[1],i,genestructure_info_dict[i]['gene_name'][0])
		print >>INTRON,"</zooms>"
	print '[genestructure_print: genestructure_info_dict detailed]',gj.print_dict(genestructure_info_dict)
	INTRON_TXT.close()
	gj.printFuncRun('genestructure_print')
	print
	return genestructure_str,genestructure_scale_str

def highlight_print(rows, out, ensembl_gene_abstract_type_dict=None, species='hs', mode='highlight', unit=1000000.0, duplex_genes=None, gene_link_gene_id_ls=None):
	gj.printFuncRun('highlight_print')
	rows_max_dict = {}
	legal_chr = map(str,range(1,23))+['X','Y','M','MT']+['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI']
	#legal_chr = ['chr'+i for i in legal_chr]
	for i in rows:
		print "[highlight_print]",i
		if not rows_max_dict.has_key(i[1]):
			rows_max_dict[i[1]] = i 
		else:
			if abs(int(i[7])-int(i[8])) > abs(int(rows_max_dict[i[1]][7])-int(rows_max_dict[i[1]][8])):
				rows_max_dict[i[1]] = i 
			else:
				pass
	for i in gene_link_gene_id_ls:
		if not rows_max_dict.has_key(i):
			print "[highlight_print: gene_link_gene_id_ls not found] %s"%(i)
	rows = rows_max_dict.values()
	chromosomes_str_ls = []
	chromosomes_len_dict = nested_dict()
	OUT_GENES = open(out+'.txt','w')
	OUT_GENES_ZOOM = open(out+'.zooms','w')
	gene_max_len = 0
	gene_min_len = 100000000000
	gene_max = ""
	gene_min = ""
	print >>OUT_GENES_ZOOM,"<zooms>"
	zooms_str_ls = []
	with open(out, 'w') as OUT:
		for n,row in enumerate(rows):
			#row = list(row)
			#row[12] = row[12].replace('chr','')
			if mode == "highlight":
				if ensembl_gene_abstract_type_dict.has_key(row[1]):
					print >>OUT,"%s%s %s %s highlight_type=gene,highlight_genetype=%s"%tuple([species,row[12]]+list(row[7:9])+[ensembl_gene_abstract_type_dict[row[1]]])
				else:
					print >>OUT,"%s%s %s %s highlight_type=gene,highlight_genetype=%s"%tuple([species,row[12]]+list(row[7:9])+[row[10]])
				duplex = "source" if row[1] in duplex_genes else "derived"
				print >>OUT_GENES,"%s%s %s %s %s duplex=%s"%tuple([species,row[12]]+list(row[7:9])+[row[2]]+[duplex])
				chromosome_str = "%s%s[a%s]:%s-%s"%(species,row[3],n,row[4]/unit,row[5]/unit)
				if row[3] != "MT":
					chromosomes_str_ls.append(chromosome_str)
			if mode == "linkgenes":
				if ensembl_gene_abstract_type_dict.has_key(row[1]):
					gene_type = ensembl_gene_abstract_type_dict[row[1]]
				else:
					gene_type = row[10]
				if row[12] not in legal_chr: continue
				print >>OUT,"%s%s %s %s gene_id=%s,gene_name=%s,gene_type=%s,tx_id=%s,id=%s"%tuple([species,row[12]]+list(row[7:9])+list(row[1:3])+[gene_type,row[6]]+['|'.join(list(map(str,row[1:])))])
				duplex = "source" if row[1] in duplex_genes else "derived"
				gene_name_pre = "***" if row[1] in duplex_genes else ""
				gene_name_sub = "***" if row[1] in duplex_genes else ""
				gene_name = row[2] if row[2] is not None else row[1] 
				if row[1] in duplex_genes and gene_name.startswith('chr'):
					gene_name_pre,gene_name_sub = '*','*'
				print >>OUT_GENES,"%s%s %s %s %s%s%s duplex=%s,gene_id=%s,gene_name=%s,gene_type=%s,tx_id=%s"%tuple([species,row[12]]+list(row[7:9])+[gene_name_pre,gene_name,gene_name_sub]+[duplex]+list(row[1:3])+[gene_type,row[6]])
				zooms_str = "<zoom> \n chr = %s \n start = %sb \n end =%sb \n scale=1 \n </zoom>"%tuple(list(row[3:6]))
				zooms_str_ls.append(zooms_str)
				if abs(int(row[7])-int(row[8])) > gene_max_len:
					gene_max_len = abs(int(row[7])-int(row[8]))
					gene_max = row
				if abs(int(row[7])-int(row[8])) < gene_min_len:
					gene_min_len = abs(int(row[7])-int(row[8]))
					gene_min = row
				print >>OUT_GENES_ZOOM,zooms_str
				chromosome_str = "%s%s[a%s]:%s-%s"%(species,row[12],n,row[7]/unit,row[8]/unit)
				chromosomes_len_dict['a%s'%(n)]['len'] = abs(int(row[7])-int(row[8]))
				chromosomes_len_dict['a%s'%(n)]['log10len'] = np.log2(abs(int(row[7])-int(row[8])))

				#if "MT" not in row[12] and "M" not in row[12]:
				chromosomes_str_ls.append(chromosome_str)
				print "[highlight_print: append to chromosomes_str_ls]",chromosome_str,row[1]
				#if "MT" in row[12] or "M" in row[12]:
				#	print "[highlight_print: chromosomes_str_ls found MT, DROP]",row
					
	OUT_GENES.close()
	print >>OUT_GENES_ZOOM,"</zooms>"
	OUT_GENES_ZOOM.close()
	print "[highlight_print: max gene len] %s"%(gene_max_len), gene_max
	print "[highlight_print: min gene len] %s"%(gene_min_len), gene_min
	s2 = 1
	l2 = gene_min_len
	with open(out+'.zooms2','w') as OUT_GENES_ZOOM:
		print >>OUT_GENES_ZOOM,"<zooms>"
		for row in rows:
			l1 = abs(int(row[7])-int(row[8]))
			#s1 = (l2*s2)/(l1*np.log10(l2-l1)) if l2 > l1 else 1 
			s1 = (l2*s2)/(l1*np.log10(l1-l2)) if l2 < l1 else 1
			print >>OUT_GENES_ZOOM,"<zoom> \n chr = %s%s \n start = %sb \n end =%sb \n scale=%s \n </zoom>"%tuple([species,row[12]]+list(row[7:9])+[s1])
			print "[highlight_print: log2 scaled coordinate based on min len]",row,l2,s2,l1,s1
		print >>OUT_GENES_ZOOM,"</zooms>"
	#if mode == "linkgenes":
	print "[highlight_print: chromosome_str]",';'.join(sorted(chromosomes_str_ls))
	print chromosomes_len_dict
	log10len_all = sum([j['log10len'] for i,j in chromosomes_len_dict.items()])
	for i,j in chromosomes_len_dict.items():
		chromosomes_len_dict[i]['degree'] = j['log10len']/log10len_all
		if j['log10len']/log10len_all == 1: 
			chromosomes_len_dict[i]['degree'] = 0.99
	print "[highlight_print: chromosomes_len_dict]",gj.print_dict(chromosomes_len_dict)
	chromosomes_scale_str = ','.join([i+':'+str(j['degree'])+'rn' for i,j in chromosomes_len_dict.items()])
	print "[highlight_print: chromosomes_scale_str(degree)]",chromosomes_scale_str
	print "[highlight_print: chromosomes_scale_str degrees]",[j['degree'] for i,j in chromosomes_len_dict.items()]
	gj.printFuncRun('highlight_print')
	print

	return ';'.join(chromosomes_str_ls),chromosomes_scale_str


def create_circos_conf(save_conf='/Users/gongjing/kuaipan/10-bioinformatics/circos-tutorials-0.67/tutorials/5/test/circos_rise.conf', species='hs', chromosome_str=None, chromosomes_scale_str=None, zooms_conf=None, highlight=None, genestructure=None, anno=None, genes=None, link=None, scale=1):

	karyotype_dict = {'hs':'karyotype.human.hg38.txt','mm':'karyotype.mouse.mm10.txt','yeast':'karyotype.yeast.chrAdd.txt'}
	chromosome_str_default_dict = { 'hs':';'.join(['hs'+str(i) for i in range(1,23)+['X','Y']]),
	                                'mm':';'.join(['mm'+str(i) for i in range(1,20)+['X','Y']]),}
	#print chromosome_str_default_dict
	chromosome_str = chromosome_str_default_dict[species] if chromosome_str is None else chromosome_str
	chromosomes_scale_str = "" if chromosomes_scale_str is None else chromosomes_scale_str

	if scale:
		link_thickness = 10
		label_size = '35p'
	else:
		link_thickness = 1
		label_size = '15p'


	gj.printFuncRun("create_circos_conf: %s"%(save_conf))

	conf_str = """<<include colors_fonts_patterns.conf>>


<<include ideogram.conf>>
<<include ticks.conf>>


<image>
<<include etc/image.conf>>
radius* = 1500p  # default=1500p
image_map_use      = yes
image_map_name     = circosmap
</image>


karyotype   = /Share/home/zhangqf5/gongjing/software/circos-0.69-4/data/karyotype/%s

chromosomes_units = 1000000
chromosomes       = %s
chromosomes_display_default = no


#chromosomes = hs12[a]:0-90;hs12[b]:90-100;hs12[c]:100-)

#chromosomes_scale = b:10
#chromosomes_scale = a18:0.25r,a6:0.25r
#chromosomes_scale = /./=1rn
chromosomes_scale = %s
chromosomes_color = /./=white_a5


# zooms
#<<include %s>>

my_color_red = 202,75,78
my_color_green = 83,169,102
my_color_blue = 74,113,178
my_color_purple = 129,112,182
my_color_brown = 205,185,111
my_color_cyan = 98,180,208
my_color_cyan_original = 0,191,255

# highlight
<highlights>
layers_overflow=collapse\
z = 100

<highlight>
url              = type=genes,id=[id]
file = %s
r0 = 0.99r
r1 = 1.01r

fill_color = conf(., my_color_cyan_original)
#my_color_red = 202,75,78


<rules>
<rule>
condition = var(gene_type) eq "protein_coding" 
fill_color = conf(., my_color_red)
</rule>
<rule>
condition = var(gene_type) eq "lncRNA"
fill_color = conf(., my_color_green)
</rule>
<rule>
condition = var(gene_type) eq "ncRNA"
fill_color = conf(., my_color_blue)
</rule>
<rule>
condition = var(gene_type) eq "snoRNA"
fill_color = conf(., my_color_purple)
</rule>
<rule>
condition = var(gene_type) eq "miRNA"
fill_color = conf(., my_color_brown)
</rule>
<rule>
condition = var(gene_type) eq "rRNA"
fill_color = conf(my_color_cyan)
</rule>
<rule>
condition = var(gene_type) eq "snRNA"
fill_color = red
</rule>
<rule>
condition = var(gene_type) eq "tRNA"
fill_color = green
</rule>
<rule>
condition = var(gene_type) eq "others"
fill_color = grey
</rule>
<rule>
condition = var(gene_type) eq "pseudogene"
fill_color = black
</rule>
<rule>
condition = var(gene_type) eq "repeats" or var(gene_type) eq "Transposable_element" or var(gene_type) eq "transposable_element" or var(gene_type) eq "Low_complexity"
fill_color = conf(., my_color_cyan_original)
</rule>
</rules>

</highlight>

<highlight>
file = %s
r0 = 0.945r
r1 = 0.955r
z = 0
fill_color = grey
</highlight>

<highlight>
url              = type=gene_structure,id=[id]
file = %s
r0 = 0.93r
#r0 = 0.99r
r1 = 0.98r
#r1 = 1.01r
fill_color = black
z = 2

<rules>
<rule>
ues = no
condition = var(sequence_type) eq "CDS"
fill_color = conf(., my_color_red) # red
#fill_color = black
</rule>
<rule>
ues = no
condition = var(sequence_type) eq "5UTR"
fill_color = conf(., my_color_blue) # blue
#fill_color = dblue
r0 = 0.94r
r1 = 0.96r
z = 5
</rule>
<rule>
ues = no
condition = var(sequence_type) eq "3UTR"
fill_color = conf(., my_color_green) # green
#fill_color = dblue
r0 = 0.94r
r1 = 0.96r
z = 5
</rule>
<rule>
ues = no
condition = var(sequence_type) eq "exon"
fill_color = lgrey
</rule>
<rule>
ues = no
condition = var(sequence_type) eq "Exon"
fill_color = lgrey
</rule>
<rule>
ues = no
condition = var(gene_type) eq "repeats" or var(sequence_type) eq "transposable_element"
fill_color = conf(., my_color_cyan_original)
</rule>
</rules>

</highlight>

</highlights>


<plots>

type            = tile
#layers_overflow = hide

# RBP binding site
<plot>

url              = type=RBP,id=[id]

<axes>
# Show axes only on ideograms that have data for this track: show = data
show = all

thickness = 3
color     = conf(., my_color_brown)
#fill_color = vorange

<axis>
spacing   = 1.1r
</axis>
</axes>

file = %s
color = conf(., my_color_brown)
fill_color = conf(., my_color_brown)
max = 1
min = 0
type  = histogram
r0    = 0.81r
r1    = 0.91r
extend_bin = no

<rules>
<rule>
use = yes
condition = var(anno_category) ne "RBP"
show = no
</rule>
</rules>

</plot>


# SNV lable
<plot>
url              = type=SNV,id=[id]
<axes>
# Show axes only on ideograms that have data for this track
show = all

thickness = 3
color     = conf(., my_color_blue)

<axis>
spacing   = 1.1r
</axis>
</axes>

file = %s
color = conf(., my_color_blue)
#fill_color = conf(., my_color_blue)
max = 1
min = 0
type  = histogram
r0    = 0.68r
r1    = 0.78r
extend_bin = no

<rules>
<rule>
use = yes
condition = var(anno_category) ne "SNV"
show = no
</rule>

<rule>
use = yes
condition = var(anno_source) eq "PanCancer"
color = conf(., my_color_blue)
fill_color = conf(., my_color_blue)
</rule>

<rule>
use = yes
condition = var(anno_source) eq "dbSNP"
color = conf(., my_color_green)
fill_color = conf(., my_color_green)
</rule>

</rules>


</plot>

# RNA modification
<plot>
url              = type=modification,id=[id]
<axes>
# Show axes only on ideograms that have data for this track
show = all

thickness = 3
color     = conf(., my_color_green) # green 

<axis>
spacing   = 1.1r
</axis>
</axes>

file = %s
color = 83,169,102
#fill_color = dgreen
max = 1
min = 0
type  = histogram
r0    = 0.56r
r1    = 0.66r
extend_bin = no

<rules>
<rule>
use = yes
condition = var(anno_category) ne "modification"
show = no
</rule>

<rule>
use = yes
condition = var(anno_source) eq "RADAR"
color = conf(., my_color_blue) # blue
fill_color = conf(., my_color_blue)
</rule>

<rule>
use = yes
condition = var(anno_source) eq "DARNED"
color = conf(., my_color_blue) # blue
fill_color = conf(., my_color_blue)
</rule>

<rule>
use = yes
condition = var(anno_type) eq "m6A"
color = conf(., my_color_red) # red
fill_color = conf(., my_color_red)
</rule>

<rule>
use = yes
condition = var(anno_type) eq "m5C" or var(anno_type) eq "N"
color = conf(., my_color_green) # green
fill_color = conf(., my_color_green)
</rule>

<rule>
use = yes
condition = var(anno_type) eq "Y"
color = conf(., my_color_purple) # purple
fill_color = conf(., my_color_purple)
</rule>


</rules>

</plot>


<plot>
type             = text
color            = dred
file             = %s

r0 = 1r
r1 = 1r+800p

show_links     = yes
link_dims      = 4p,4p,8p,4p,4p
link_thickness = 2p
link_color     = red

label_size   = %s
label_font   = condensed
label_rotate = yes

padding  = 0p
rpadding = 0p

label_snuggle             = yes  # snuggle: settle or move into a warm, comfortable position.
#snuggle_refine                 = yes
#snuggle_tolerance               = 0.3r


<rules>
<rule>
use = yes
condition = var(duplex) eq "source"
color = black
label_font = bold
label_size = 45p
</rule>

<rule>
condition = var(gene_type) eq "protein_coding" 
color = conf(., my_color_red)
</rule>
<rule>
condition = var(gene_type) eq "lncRNA"
color = conf(., my_color_green)
</rule>
<rule>
condition = var(gene_type) eq "ncRNA"
color = conf(., my_color_blue)
</rule>
<rule>
condition = var(gene_type) eq "snoRNA"
color = conf(., my_color_purple)
</rule>
<rule>
condition = var(gene_type) eq "miRNA"
color = conf(., my_color_brown)
</rule>
<rule>
condition = var(gene_type) eq "rRNA"
color = conf(., my_color_cyan)
</rule>
<rule>
condition = var(gene_type) eq "snRNA"
color = red
</rule>
<rule>
condition = var(gene_type) eq "tRNA"
color = green
</rule>
<rule>
condition = var(gene_type) eq "others"
color = grey
</rule>
<rule>
condition = var(gene_type) eq "pseudogene"
color = black
</rule>
<rule>
condition = var(gene_type) eq "repeats" or var(gene_type) eq "Transposable_element"
color = conf(., my_color_cyan_original)
</rule>

</rules>

</plot>


</plots>

<links>

radius = 0.54r
crest  = 1
ribbon = yes
color  = grey

bezier_radius        = 0r
#bezier_radius_purity = 0.5
stroke_color = grey
stroke_thickness = %s
thickness    = 10  
category = all
method = all
duplex = all


<link>
#file         = data/5/RISE0256681.link
url              = type=interaction,id=[id]
file = %s

<rules>

<rule>
use = yes
condition = var(duplex) eq "source"
thickness = 15
color = black
stroke_color = black
z=10  # control relative layer
</rule>

<rule>
use = yes
condition = var(category) eq "GlobalAnalysis"
color = lred
stroke_color = lred
#condition = var(method) eq "PARIS"
</rule>

<rule>
use = yes
condition = var(category) eq "TargetedAnalysis"
color = dblue
stroke_color = dblue
#thickness = 8
z=6
</rule>

<rule>
use = yes
condition = var(category) eq "database"
color = dgreen
stroke_color = dgreen
#thickness = 8
z=8
</rule>


</rules>

</link>

</links>

<<include etc/housekeeping.conf>>
data_out_of_range* = trim

	"""%(karyotype_dict[species], chromosome_str, chromosomes_scale_str, zooms_conf, highlight, highlight, genestructure, anno, anno, anno, genes, label_size, link_thickness, link)

	with open(save_conf,'w') as CONF:
		print >>CONF,conf_str
	gj.printFuncRun("create_circos_conf: finished")
	print 

def id_all_track_coordinate_convert(id='ENSG00000100320', dir=None,chromosomes_units=1000000.0, species='hs'):
	if dir is None:
		dir = '/Share/home/zhangqf5/gongjing/DNA-RNA-Protein-interaction-correlation-12-18/results/overlap/duplex_test'

	linkgenes = dir+'/'+id+'.linkgenes'
	linkgenes_txt = dir+'/'+id+'.linkgenes.txt'
	link = dir+'/'+id+'.link'
	genestructure = dir+'/'+id+'.genestructure'
	genestructure_intron = dir+'/'+id+'.genestructure.introns'
	genestructure_intron_txt = dir+'/'+id+'.genestructure.introns.txt'
	anno = dir+'/'+id+'.anno'

	gj.printFuncRun('id_all_track_coordinate_convert')
	gj.printFuncArgs()

	linkgenes_dict = read_linkgenes(linkgenes)
	genestructure_dict,generegion_dict = read_genestructure(genestructure_intron_txt)

	convert_ls = [linkgenes, linkgenes_txt, link, genestructure, anno]
	idx_dict = {0:'linkgenes',1:'linkgenes_txt',2:'link',3:'genestructure',4:'anno'}
	txt_info_dict = nested_dict(1, list)
	chromosomes_str_ls = []
	chromosomes_len_dict = nested_dict()
	sequence_type_color_dict = {'CDS':'#C44E52','5UTR':'#4C72B0','3UTR':'#55A868','exon':'#A9A9A9','Exon':'#A9A9A9','Transposable_element':'#00BFFF','transposable_element':'#00BFFF','repeats':'#00BFFF'}
	"""
	anno_source_color_dict = {'human_RADAR':'#4C72B0','human_DARNED':'#4C72B0','mouse_RADAR':'#4C72B0','mouse_DARNED':'#4C72B0',
                              'DARNED':'#4C72B0','RADAR':'#4C72B0',
                              'human_editing':'#4C72B0','mouse_editing':'#4C72B0',
							  'human_RMBasem6A':'#C44E52','human_RMBasem5C':'#55A868','human_RMBasePseudoU':'#8172B2','human_RMBaseOthers':'#CCB974','human_RMBaseNm':'#64B5CD',
							  'mouse_RMBasem6A':'#C44E52','mouse_RMBasem5C':'#55A868','mouse_RMBasePseudoU':'#8172B2',
                              'human_modification':'#C44E52','mouse_modification':'#C44E52',
							  'human_pancan':'#4C72B0','human_ucscSnp142':'#55A868','mouse_ucscSnp142':'#55A868',
							  'human_clipdb':'#CCB974','mouse_clipdb':'#CCB974'}
	"""
	anno_source_color_dict = {'CLIPDB':{'RBP':'#CCB974'},
                              'DARNED':{'editing':'#4C72B0'},
                              'RADAR':{'editing':'#4C72B0'},
                              'dbSNP':{'SNP':'#55A868'},
                              'PanCancer':{'mutation':'#4C72B0'},
                              'RMBase':{'m6A':'#C44E52', 
                                        'm5C':'#55A868','N':'#55A868',
                                        'Y':'#8172B2', 
                                        'Am':'#64B5CD','Cm':'#64B5CD','Gm':'#64B5CD','m5Um':'#64B5CD','Um':'#64B5CD', 
                                        }}
	modOther_ls = ['ac4c','acp3U','D','f5C','galQtRNA','I','i6A','m1A','m1G','m1I','m1Y','m2,2G','m2_2G','m2G','m3C','m5U','m7G','o2yW','QtRNA','t6A','tm5s2U','tm5U','xG','xU','Ym','yW']
	for mod in modOther_ls:
		anno_source_color_dict['RMBase'][mod] = '#CCB974'

	BioCircosGenome_str_ls = []
	for n,txt in enumerate(convert_ls):
		txt_type = idx_dict[n]
		TXT_NEW = open(txt+'.convert','w')
		with open(txt,'r') as TXT:
			gj.printFuncRun("id_all_track_coordinate_convert: process => file %s: %s"%(n,txt))
			for n_line,line in enumerate(TXT):
				line = line.strip()
				#print " - txt line %s:"%(n_line),line
				
				if not line or line.startswith('#'): continue
				arr = line.split(' ')
				info_dict = {i.split('=')[0]:i.split('=')[1] for i in arr[-1].split(',')}
				if n == 2:
					chr1,start1,end1,chr2,start2,end2 = arr[0:6]
					start1,end1 = sorted(map(int,[start1,end1]))
					start2,end2 = sorted(map(int,[start2,end2]))
					gene_id1 = info_dict['gene_id1']
					gene_id2 = info_dict['gene_id2']
					if not genestructure_dict.has_key(gene_id1) or not genestructure_dict.has_key(gene_id2): continue
					gene_name1 = genestructure_dict[gene_id1]['gene_name'][0]
					gene_name2 = genestructure_dict[gene_id2]['gene_name'][0]
					#gene_name1 = info_dict['gene_name1']
					#gene_name2 = info_dict['gene_name2']
					start1_new,end1_new,region1_convert_code = gene_id_region_region_convert(genestructure_dict=genestructure_dict, gene_id=gene_id1, start=start1, end=end1)
					start2_new,end2_new,region2_convert_code = gene_id_region_region_convert(genestructure_dict=genestructure_dict, gene_id=gene_id2, start=start2, end=end2)
					print >>TXT_NEW,' '.join(map(str,[chr1,start1_new,end1_new,chr2,start2_new,end2_new])+arr[6:])
					print "- txt line %s:"%(n_line),line,start1_new,end1_new,region1_convert_code,start2_new,end2_new,region2_convert_code
					#print >>TXT_NEW,' '.join(map(str,[gene_name1,start1_new,end1_new,gene_name2,start2_new,end2_new])+arr[6:])
					biocircos_link_dict = OrderedDict()
					biocircos_link_dict['fusion'] = "%s (%s)"%(info_dict['category'], info_dict['method'])
					biocircos_link_dict['g1chr'] = "%s"%(gene_name1)
					biocircos_link_dict['g1start'] = "%s"%(start1_new)
					biocircos_link_dict['g1end'] = "%s"%(end1_new)
					biocircos_link_dict['g1name'] = "%s"%(gene_name1)
					biocircos_link_dict['g2chr'] = "%s"%(gene_name2)
					biocircos_link_dict['g2start'] = "%s"%(start2_new)
					biocircos_link_dict['g2end'] = "%s"%(end2_new)
					biocircos_link_dict['g2name'] = "%s"%(gene_name2)
					#biocircos_link_dict = {"fusion": "%s--%s"%(gene_name1,gene_name2), "g1chr": "%s"%(gene_name1), "g1start": "%s"%(start1_new), "g1end": "%s"%(end1_new), "g1name": "%s"%(gene_name1), "g2chr": "%s"%(gene_name2), "g2start": "%s"%(start2_new), "g2end": "%s"%(end2_new), "g2name": "%s"%(gene_name2)}
					txt_info_dict[txt_type].append(biocircos_link_dict)
					txt_info_dict['link_category_type'].append(info_dict['category'])
					"""
					gene_id1_start = linkgenes_dict[gene_id1]['start']
					gene_id2_start = linkgenes_dict[gene_id2]['start']
					start1 = gene_id1_start + np.log2(start1 - gene_id1_start) if start1 > gene_id1_start else gene_id1_start
					end1 = start1 + np.log2(abs(int(arr[2])-int(arr[1])))
					start2 = gene_id2_start + np.log2(start2 - gene_id2_start) if start2 > gene_id2_start else gene_id2_start
					end2 = start2 + np.log2(abs(int(arr[5])-int(arr[4])))
					print >>TXT_NEW,' '.join(map(str,[chr1,int(ceil(start1)),int(ceil(end1)),chr2,int(ceil(start2)),int(ceil(end2))])+arr[6:])
					"""
				else:
					chr1,start1,end1 = arr[0:3]
					start1,end1 = sorted(map(int,[start1,end1]))
					if n == 4:
						gene_id = info_dict['anno_gene_info'].split('|')[6]
						if not genestructure_dict.has_key(gene_id): continue
						gene_name = genestructure_dict[gene_id]['gene_name'][0]
						#gene_name = info_dict['anno_gene_info'].split('|')[7]
					else:
						gene_id = info_dict['gene_id']
						if not genestructure_dict.has_key(gene_id): continue
						gene_name = genestructure_dict[gene_id]['gene_name'][0]
						#gene_name = info_dict['gene_name']
					if n == 1 or n == 0:
						start1_new = int(genestructure_dict[gene_id]['start_log2'])
						end1_new = int(genestructure_dict[gene_id]['end_log2'])
						region1_convert_code = "assign_directly"
					else:
						start1_new,end1_new,region1_convert_code = gene_id_region_region_convert(genestructure_dict=genestructure_dict, gene_id=gene_id, start=start1, end=end1)
					print >>TXT_NEW,' '.join(map(str,[chr1,start1_new,end1_new])+arr[3:])
					print "- txt line %s:"%(n_line),line,start1_new,end1_new,region1_convert_code
					#print >>TXT_NEW,' '.join(map(str,[gene_name,start1_new,end1_new])+arr[3:])
					"""
					gene_id_start = linkgenes_dict[gene_id]['start']
					start1 = gene_id_start + np.log2(start1 - gene_id_start) if start1 > gene_id_start else gene_id_start
					#end1 = gene_id_start + np.log2(end1 - gene_id_start) if end1 > gene_id_start else gene_id_start
					end1 = start1 + np.log2(abs(int(arr[2])-int(arr[1])))
					print >>TXT_NEW,' '.join(map(str,[chr1,int(ceil(start1)),int(ceil(end1))])+arr[3:])
					"""
					if n == 1:
						#if "MT" not in chr1:
						chromosomes_str_ls.append("%s[a%s]:%.6f-%.6f"%(chr1,n_line,start1_new/chromosomes_units,end1_new/chromosomes_units))
						chromosomes_len_dict['a%s'%(n_line)]['len'] = abs(start1_new-end1_new)
						chromosomes_len_dict['a%s'%(n_line)]['log10len'] = np.log2(abs(start1_new-end1_new)) if start1_new != end1_new else 1
					biocircos_link_dict = OrderedDict()
					if n == 3:
						biocircos_link_dict['chr'] = "%s"%(gene_name)
						biocircos_link_dict['start'] = "%s"%(start1_new)
						biocircos_link_dict['end'] = "%s"%(end1_new)
						biocircos_link_dict['color'] = "%s"%(sequence_type_color_dict[info_dict['sequence_type']])
						biocircos_link_dict['des'] = "%s"%(info_dict['sequence_type'])
						txt_info_dict[txt_type].append(biocircos_link_dict)
						txt_info_dict['genestructure_type'].append(info_dict['sequence_type'])
					if n == 4:
						biocircos_link_dict['chr'] = "%s"%(gene_name)
						biocircos_link_dict['start'] = "%s"%(start1_new)
						biocircos_link_dict['end'] = "%s"%(end1_new)
						biocircos_link_dict['color'] = "%s"%(anno_source_color_dict[info_dict['anno_source']][info_dict['anno_type']])
						biocircos_link_dict['des'] = "%s;%s"%(info_dict['anno_type'],info_dict['anno_source'])
						#biocircos_link_dict['name'] = "%s"%(info_dict['anno_source'])
						#biocircos_link_dict['value'] = "%s"%(info_dict['anno_type'])
						txt_info_dict[txt_type].append(biocircos_link_dict)
						txt_info_dict['anno_type'].append(info_dict['anno_type'])
					if n == 0:
						BioCircosGenome_str_ls.append([gene_name, end1_new])

		TXT_NEW.close()

	save_conf = '/Share/home/zhangqf5/gongjing/software/circos-tutorials-0.67/tutorials/5/test/rise_conf/%s.convert.conf'%(id)
	chromosomes_str = ';'.join(chromosomes_str_ls)
	print "[id_all_track_coordinate_convert: chromosomes_str of circos]",chromosomes_str
	log10len_all = sum([j['log10len'] for i,j in chromosomes_len_dict.items()])
	for i,j in chromosomes_len_dict.items():
		chromosomes_len_dict[i]['degree'] = j['log10len']/log10len_all
		if j['log10len']/log10len_all == 1: chromosomes_len_dict[i]['degree'] = 0.99
	print "[id_all_track_coordinate_convert: chromosomes_len_dict of circos]",chromosomes_len_dict
	chromosomes_scale_str = ','.join([i+':'+str(j['degree'])+'rn' for i,j in chromosomes_len_dict.items()])
	print "[id_all_track_coordinate_convert: chromosomes_scale_str of circos]",chromosomes_scale_str

	create_circos_conf(save_conf=save_conf, species=species, chromosome_str=chromosomes_str, chromosomes_scale_str=chromosomes_scale_str, zooms_conf=None, highlight=linkgenes+'.convert', genestructure=genestructure+'.convert', anno=anno+'.convert', genes=linkgenes_txt+'.convert', link=link+'.convert')
	circos_plot(conf=save_conf,savefn_dir=dir,savefn=id+'.convert.png')
	gj.printFuncRun('id_all_track_coordinate_convert')



def gene_id_region_region_convert(genestructure_dict=None, gene_id=None, start=None, end=None):
	#if gene_id == "ENSG00000197852":
		#print gene_id, start, end 
	start,end = sorted([start, end])
	if start == ".":  # deal with entry like: hs1 71067716 71067695 hs6 26157105 .
		start = int(end)-1
	region_len = abs(int(end) - int(start))
	status = 0
	convert_code = 0
	for i,j in zip(genestructure_dict[gene_id]['regions'],genestructure_dict[gene_id]['regions_log2']):
		#if gene_id == "ENSG00000197852": print i,j
		if start == genestructure_dict[gene_id]['start'] and end == genestructure_dict[gene_id]['end']: # for gene ha only one region: from start to end
			start, end = genestructure_dict[gene_id]['start_log2'],genestructure_dict[gene_id]['end_log2']
			status += 1
			convert_code = 1
			#if gene_id == "ENSG00000197852": print 'converted',start,end
			break
		if start == i[0] and end == i[1]: # start,end exactly same as region
			start = j[0]
			end = j[1]
			#if gene_id == "ENSG00000197852": print 'converted',start,end
			status = 1
			convert_code = 2
			break
		elif (start > i[0] and end < i[1]) or (start > i[0] and end == i[1]) or (start == i[0] and end < i[1]): # start,end as part of a gene's region
			start = j[0] if int(start) == i[0] else j[0] + np.log2(abs(int(start)-i[0])) # if start=i[0] => log0 => -inf
			#end = start + np.log2(region_len)
			end = j[0] + np.log2(abs(int(end))-i[0])
			status = 1
			convert_code = 3
			#if gene_id == "ENSG00000197852": print 'converted',start,end
			break
		elif start < i[0] and end < i[0]: # start,end all small than a gene's region, assign the 1st gene's region start,end
			start = j[0]
			end = j[1]
			status += 1
			convert_code = 4
			break
		elif start < i[0] and end > i[0]: # start out of region, end within region, set start as region start
			start = j[0]
			end = j[0] + np.log2(end-i[0])
			status += 1
			convert_code = 5
			break
		else:
			pass
			#if gene_id == 'ENSG00000204345':
	if status == 0: # if start,end out of gene's all regions(too small/big), set as the gene's region middle point
			print "[error]non process:",gene_id,start,end,genestructure_dict[gene_id]
			start = int((genestructure_dict[gene_id]['start_log2']+genestructure_dict[gene_id]['end_log2'])/2.0)
			end = start 
			convert_code = 6
			#sys.exit()
	#if int(ceil(start)) == int(ceil(end)):
	#	return int(ceil(start)),int(ceil(end))+1
	if ceil(end) > genestructure_dict[gene_id]['end_log2'] or ceil(start) < genestructure_dict[gene_id]['start_log2']:
		print "[error]process out of range:",gene_id,start,end,genestructure_dict[gene_id]
		convert_code = 7
		#sys.exit()
	if int(ceil(start)) == genestructure_dict[gene_id]['end_log2']:
		start, end = int(ceil(start))-1,int(ceil(end)) 
		convert_code = 8
	#if start == end:
		#end = end + 1
	return int(ceil(start)),int(ceil(end)),convert_code

def read_linkgenes(linkgenes):
	gj.printFuncRun('read_linkgenes')
	gj.printFuncArgs()
	linkgenes_dict = nested_dict()
	with open(linkgenes,'r') as LINKGENES:
		for line in LINKGENES:
			if not line or line.startswith('#'): continue
			chromosome,start,end,info = line.strip().split(' ')
			info_dict = {i.split('=')[0]:i.split('=')[1] for i in info.split(',')}
			gene_id = info_dict['gene_id']
			for i,j in info_dict.items():
				linkgenes_dict[gene_id][i] = j
			linkgenes_dict[gene_id]['chr'] = chromosome
			linkgenes_dict[gene_id]['start'] = int(start)
			linkgenes_dict[gene_id]['end'] = int(end)
	print "[read_linkgenes: linkgenes_dict]",linkgenes_dict
	print 
	return linkgenes_dict

def read_genestructure(genestructure_txt='/Users/gongjing/kuaipan/1-labs-Zhangqf/10-DNA-RNA-protein-interaction-corrrelation/results/overlap/duplex/ENSG00000100320.genestructure.introns.txt'):
	gj.printFuncRun('read_genestructure')
	gj.printFuncArgs()
	genestructure_dict = nested_dict(2,list)
	generegion_dict = nested_dict()
	with open(genestructure_txt,'r') as TXT:
		for line in TXT:
			line = line.strip()
			if not line or line.startswith('#'): continue
			arr = line.split(' ')
			info_dict = {i.split('=')[0]:i.split('=')[1] for i in arr[-1].split(',')}
			gene_id = info_dict['gene_id']
			genestructure_dict[gene_id]['gene_id'] = [gene_id] 
			genestructure_dict[gene_id]['gene_name'] = [info_dict['gene_name']] 
			generegion_dict[arr[1]+'_'+arr[2]]['type'] = info_dict['type']
			genestructure_dict[gene_id]['regions'].append(sorted([int(arr[1]),int(arr[2])]))
	for i,j in genestructure_dict.items():
		genestructure_dict[i]['regions'] = sorted(j['regions'])
		genestructure_dict[i]['start'] = min([q for p in j['regions'] for q in p])
		genestructure_dict[i]['end'] = max([q for p in j['regions'] for q in p])
		regions_log2 = []
		for n,k in enumerate(sorted(j['regions'])):
			if n == 0:
				start = float(k[0])  # set as tx start
				#start = 1 # set as 1
				end = start + np.log2(abs(int(k[0])-int(k[1]))) if k[0] != k[1] else start
			else:
				start = regions_log2[n-1][1]
				end = start + np.log2(abs(int(k[0])-int(k[1]))) if k[0] != k[1] else start
			regions_log2.append([ceil(start),ceil(end)])
		genestructure_dict[i]['regions_log2'] = regions_log2
		genestructure_dict[i]['start_log2'] = min([q for p in j['regions_log2'] for q in p])
		genestructure_dict[i]['end_log2'] = max([q for p in j['regions_log2'] for q in p])
	print "[read_genestructure: genestructure_dict]",genestructure_dict
	#gj.print_dict(generegion_dict)
	#print genestructure_dict['ENSG00000153037']
	gj.printFuncRun('read_genestructure')
	print 
	return genestructure_dict,generegion_dict

def parse_html(html='/Share/home/zhangqf5/gongjing/DNA-RNA-Protein-interaction-correlation-12-18/results/overlap/duplex/ENSG00000100320.convert.html', savefn_html=None):
	gj.printFuncRun('parse_html')
	gj.printFuncArgs()
	if savefn_html is None:
		savefn = html.replace('html', 'resize.html')
	else:
		savefn = savefn_html
	HTML_OUT = open(savefn, 'w')
	html_prefix = """<html>
<head>
	<meta charset="utf-8">
	<title>mycircos</title>
	<style>
	* {
  	padding: 0;
  	margin: 0;
	}
	.fit{%s}
	</style>
	<script type="text/javascript" src="https://ajax.googleapis.com/ajax/libs/jquery/1.9.0/jquery.min.js"></script>
	<script type="text/javascript" src="http://davidlynch.org/projects/maphilight/jquery.maphilight.min.js"></script>
	<script type="text/javascript" src="http://davidjbradshaw.com/imagemap-resizer/js/imageMapResizer.min.js"></script>
	<script type="text/javascript">$(function() {
		//$('map').imageMapResize();
		$('.map').maphilight();
		//$('area').mouseenter(function() {
  		//	$( "area" ).append( "<div>Handler for .mouseenter() called.</div>" );	
		//});
		$('area').mouseleave(function(d) {
  			
  			$("#BioCircosLINKTooltip").css("display", "none");

  			 
		});
		$('area').mouseenter(function(d) {
  			//$( "area" ).append( "<div>Handler for .mouseover() called."+$( this ).attr("title")+"x:"+d.pageX+"</div>" );
  			$("#BioCircosLINKTooltip").css("display", "block");
  			$("#BioCircosLINKTooltip").css("left", (d.pageX) + "px")
                       .css("top", (d.pageY + 20) + "px");
            var content = content_format($( this ).attr("alt"))
            $("#BioCircosLINKTooltip").html(content);

  			 
		});
	});</script>

<script type="text/javascript">
function content_format(str) {
	var type = str.match(/type=.*,/);
    var anno_source = str.match(/anno_source=.*,/);
	var id = str.match(/id=.*/);
	var id_arr = id[0].split('=')[1].split('|');

	var interaction_ls = ["chr1", "start1", "end1", "chr2", "start2", "end2", "ID", "score", "strand1", "strand2", "ensembl_gene_id1", "ensembl_gene_name1", "ensembl_gene_id2", "ensembl_gene_name2", "ensembl_gene_type1", "ensembl_gene_type2", "ensembl_gene_abstract_type1","ensembl_gene_abstract_type2","method","species","cell_line"] ;
	var gene_ls = ["gene_id","gene_name","gene_start","gene_end","gene_type","tx_id","tx_start","tx_end","tx_type","tx_abstract_type","feature_type","chr","feature_start","feature_end"]
	var gene_structure_ls = ["gene_id","gene_name","gene_start","gene_end","gene_type","tx_id","tx_start","tx_end","tx_type","tx_abstract_type","feature_type","chr","feature_start","feature_end"]
	//var modification_ls = ["chr1","start1","end1","RISE_id","score1","strand1","gene_id1","gene_name1","gene_type1","gene_abstract_type1","chr", "start", "end", "id", "score", "strand"]
	var RMBase_ls = ['chr1', 'start1', 'end1', 'RISE_id', 'score1', 'strand1', 'gene_id1', 'gene_name1', 'gene_type1', 'gene_abstract_type1', 'chr', 'start', 'end', 'id', 'score', 'strand', 'id2', 'type', 'supportNum', 'supportList', 'pubmed_id', 'gene_name', 'gene_type', 'region_type', 'sequence','category', 'source', 'reference', 'pos_gene_id', 'pos_tx_id', 'pos_genomic_context', 'overlap_base_num']
    var editing_ls = ['chr1', 'start1', 'end1', 'RISE_id', 'score1', 'strand1', 'gene_id1', 'gene_name1', 'gene_type1', 'gene_abstract_type1', 'chr', 'start', 'end', 'id', 'score', 'strand', 'category', 'source', 'type','pubmed_id','reference', 'pos_gene_id', 'pos_tx_id', 'pos_genomic_context', 'overlap_base_num']
	//var RBP_ls = ["chr1","start1","end1","RISE_id","score1","strand1","gene_id1","gene_name1","gene_type1","gene_abstract_type1",'chr', 'start', 'end', 'id', 'score', 'strand', 'gene_name', 'methods', 'cell_line', 'sampleID', 'score2']
	var RBP_ls = ['chr1', 'start1', 'end1', 'RISE_id', 'score1', 'strand1', 'gene_id1', 'gene_name1', 'gene_type1', 'gene_abstract_type1', 'chr', 'start', 'end', 'id', 'score', 'strand', 'gene_name', 'methods', 'cell_line', 'sampleID', 'score2', 'pubmed_id', 'reference', 'pos_gene_id', 'pos_tx_id', 'pos_genomic_context', 'overlap_base_num','method(experiment)','method(computation)','species']
	//var SNV_ls = ["chr1","start1","end1","RISE_id","score1","strand1","gene_id1","gene_name1","gene_type1","gene_abstract_type1",'chr', 'start', 'end', 'id', 'score', 'strand','ref_base','alt_base']
	var SNV_pancan_ls = ['chr1', 'start1', 'end1', 'RISE_id', 'score1', 'strand1', 'gene_id1', 'gene_name1', 'gene_type1', 'gene_abstract_type1', 'chr', 'start', 'end', 'id', 'score', 'strand', 'ref_base', 'alt_base', 'sample', 'cancer', 'category', 'source', 'pubmed_id', 'reference', 'pos_gene_id', 'pos_tx_id', 'pos_genomic_context', 'overlap_base_num','cancer(full)']
	var SNV_dbsnp_ls = ['chr1', 'start1', 'end1', 'RISE_id', 'score1', 'strand1', 'gene_id1', 'gene_name1', 'gene_type1', 'gene_abstract_type1', 'chr', 'start', 'end', 'id', 'score', 'strand', 'ref_base', 'alt_base', 'category', 'source', 'pubmed_id', 'reference', 'pos_gene_id', 'pos_tx_id', 'pos_genomic_context', 'overlap_base_num']
	var id_dict = {};
	var id_dict_str = "";

	if (type == "type=interaction,") {
		for (i = 0; i < id_arr.length; i += 1){id_dict[interaction_ls[i]] = id_arr[i];};
		//for (key in id_dict){id_dict_str += key+": "+id_dict[key]+"<br>"};
		id_dict_str = "Track: interaction"+'<br>'+'RNA1 name: '+id_dict['ensembl_gene_name1']+'<br>'+'RNA1 ID: '+ id_dict['ensembl_gene_id1']+'<br>'+'RNA1 gene type: '+id_dict['ensembl_gene_abstract_type1']+'<br>'+'RNA1 region: '+id_dict['chr1']+':'+id_dict['start1']+'-'+id_dict['end1']+'<br>'+'RNA2 name: '+id_dict['ensembl_gene_name2']+'<br>'+'RNA2 ID: '+ id_dict['ensembl_gene_id2']+'<br>'+'RNA2 gene type: '+id_dict['ensembl_gene_abstract_type2']+'<br>'+'RNA2 region: '+id_dict['chr2']+':'+id_dict['start2']+'-'+id_dict['end2']+'<br>'+'Method: '+id_dict['method']+'<br>'+'Species: '+id_dict['species']+'<br>'+'Cell line: '+id_dict['cell_line']
		return id_dict_str
	} else if (type == "type=genes,") {
		for (i = 0; i< id_arr.length; i += 1){id_dict[gene_ls[i]] = id_arr[i]};
		//for (key in id_dict){id_dict_str += key+": "+id_dict[key]+"<br>"};
		id_dict_str = "Track: gene"+'<br>'+'Gene name: '+id_dict['gene_name']+'<br>'+'Gene ID: '+id_dict['gene_id']+'<br>'+'Gene type: '+id_dict['tx_abstract_type']+'<br>'+'Position: chr'+id_dict['chr']+':'+id_dict['gene_start']+'-'+id_dict['gene_end']
		return id_dict_str
	} else if (type == "type=gene_structure,"){
		for (i = 0; i< id_arr.length; i += 1){id_dict[gene_structure_ls[i]] = id_arr[i]};
		//for (key in id_dict){id_dict_str += key+": "+id_dict[key]+"<br>"};
		id_dict_str = "Track: gene structure"+'<br>'+'Gene name: '+id_dict['gene_name']+'<br>'+'Sequence type: '+id_dict['feature_type']+'<br>'+'Position: chr'+id_dict['chr']+':'+id_dict['feature_start']+'-'+id_dict['feature_end']
		return id_dict_str
	} else if (type == "type=modification,"){
        if (id_arr.length == RMBase_ls.length) {
        for (i = 0; i< RMBase_ls.length; i += 1){id_dict[RMBase_ls[i]] = id_arr[i]};
        id_dict_str = "Track: processing"+'<br>'+'Category: RNA modification'+'<br>'+'Type: '+id_dict['type']+'<br>'+'Position: '+id_dict['chr']+':'+id_dict['start']+'-'+id_dict['end']+'<br>'+'Source: '+id_dict['source']
        } else {
        for (i = 0; i< editing_ls.length; i += 1){id_dict[editing_ls[i]] = id_arr[i]};
        id_dict_str = "Track: processing"+'<br>'+'Category: RNA editing'+'<br>'+'Type: '+id_dict['type']+'<br>'+'Position: '+id_dict['chr']+':'+id_dict['start']+'-'+id_dict['end']+'<br>'+'Source: '+id_dict['source']
        }
		//for (i = 0; i< modification_ls.length; i += 1){id_dict[modification_ls[i]] = id_arr[i]};
		//for (key in id_dict){id_dict_str += key+": "+id_dict[key]+"<br>"};
		//id_dict_str = "Track: modification"+'<br>'+'Type: '+id_dict['type']+'<br>'+'Position: '+id_dict['chr']+':'+id_dict['start']+'-'+id_dict['end']+'<br>'+'Source: '+id_dict['source']
		return id_dict_str
	} else if (type == "type=RBP,") {
		for (i = 0; i< RBP_ls.length; i += 1){id_dict[RBP_ls[i]] = id_arr[i]};
		//for (key in id_dict){id_dict_str += key+": "+id_dict[key]+"<br>"};
		id_dict_str = "Track: RBP binding"+'<br>'+"Protein: "+id_dict["gene_name"]+"<br>"+"Method(experiment): "+id_dict['method(experiment)']+'<br>'+"Method(computation): "+id_dict['method(computation)']+'<br>'+"Cell line: "+id_dict['cell_line']+'<br>'+'Binding score: '+parseFloat(id_dict['score2']).toFixed(3)+'<br>'+'Position: '+id_dict['chr']+':'+id_dict['start']+'-'+id_dict['end']+'<br>'+'Source: CLIPDB'
		return id_dict_str
	} else if (type == "type=SNV,") {
	    if (id_arr.length == SNV_pancan_ls.length) {
	    	for (i = 0; i< SNV_pancan_ls.length; i += 1){id_dict[SNV_pancan_ls[i]] = id_arr[i]};
	    	id_dict_str = "Track: SNV(cancer mutation)"+'<br>'+"Gene name: "+id_dict["gene_name1"]+"<br>"+"Ref base: "+id_dict["ref_base"]+"<br>"+"Alt base: "+id_dict['alt_base']+"<br>"+"ID: "+id_dict['id']+"<br>"+"Position: "+id_dict['chr']+':'+id_dict["start"]+'-'+id_dict['end']+'<br>'+'Source: Pan-Cancer'+'<br>'+'Cancer type: '+id_dict['cancer'].toUpperCase()
	    } else {
	    	for (i = 0; i< SNV_dbsnp_ls.length; i += 1){id_dict[SNV_dbsnp_ls[i]] = id_arr[i]};
	    	id_dict_str = "Track: SNV(SNP)"+'<br>'+"Gene name: "+id_dict["gene_name1"]+"<br>"+"Ref base: "+id_dict["ref_base"]+"<br>"+"Alt base: "+id_dict['alt_base']+"<br>"+"ID: "+id_dict['id']+"<br>"+"Position: "+id_dict['chr']+':'+id_dict["start"]+'-'+id_dict['end']+'<br>'+'Source: '+id_dict['source']
	    }
		//for (i = 0; i< SNV_ls.length; i += 1){id_dict[SNV_ls[i]] = id_arr[i]};
		//for (key in id_dict){id_dict_str += key+": "+id_dict[key]+"<br>"};
		//id_dict_str = "Track: SNV"+'<br>'+"Gene name: "+id_dict["gene_name1"]+"<br>"+"Ref base: "+id_dict["ref_base"]+"<br>"+"Alt base: "+id_dict['alt_base']+"<br>"+"ID: "+id_dict['id']+"<br>"+"Position: "+id_dict['chr']+':'+id_dict["start"]+'-'+id_dict['end']+'<br>'+'Source: '+id_dict['source']
		return id_dict_str
	} else {
		return "jdnqwkjdhwedh<br>ksldjqwoidje<br>klwjqdk"
	}

}
</script>

</head>

<body>

<img class="map" src="%s" usemap="#circosmap">


	"""%("width: 100%;\nheight:100%;",html.split('/')[-1].replace('html','svg'))
	print >>HTML_OUT,html_prefix



	track_base_ls = ['RBP', 'SNV', 'modification']
	track_region_ls = ['gene_structure', 'genes']
	track_resize_ls = ['interaction']
	with open(html, 'r') as HTML_IN:
		for line in HTML_IN:
			line = line.strip()
			if not line.startswith('<area'):
				print >>HTML_OUT,line.strip()
				continue
			arr = line.split(' ')
			point_num = len(arr[2].split("'")[1].split(','))/2
			#print arr 
			area_type = arr[5].split(',')[0].split('=')[-1]
			if area_type in track_base_ls and point_num == 2:
				arr[1] = arr[1].replace('poly','poly')
				rect_ls = map(int, arr[2].split("'")[1].split(','))
				rect_resize_ls = rect_resize(rect_ls, plus=1)
				arr[2] = "coords="+','.join(map(str,rect_resize_ls))+"'"
				print >>HTML_OUT,' '.join(arr[0:3]+[arr[-1]]).replace('title','alt')
			elif area_type in track_base_ls and point_num > 2:
				print >>HTML_OUT,' '.join(arr[0:3]+[arr[-1]]).replace('title','alt')
			elif area_type in track_region_ls and point_num == 2:
				arr[1] = arr[1].replace('poly','rect')
				print >>HTML_OUT,' '.join(arr[0:3]+[arr[-1]]).replace('title','alt')
			elif area_type in track_region_ls and point_num > 2:
				print >>HTML_OUT,' '.join(arr[0:3]+[arr[-1]]).replace('title','alt')
			elif area_type in track_resize_ls and point_num > 2:
				coord_ls = map(int, arr[2].split("'")[1].split(','))
				coord_resize_ls = coordinate_resize(coord_ls = coord_ls)
				arr[2] = "coords="+','.join(map(str,coord_resize_ls))+"'"
				print >>HTML_OUT,' '.join(arr[0:3]+[arr[-1]]).replace('title','alt')
			else:
				print "[unknown type or point num]: %s, %s"%(area_type, point_num), arr 
				#print >>HTML_OUT,' '.join(arr)
				sys.exit()

	html_subfix="""<div class="BioCircosLINKTooltip" id="BioCircosLINKTooltip" style="display: none; opacity: 0.8; position: absolute; background-color: white; border-style: solid; padding: 3px; border-radius: 3px;font-size:50px"></div>



</body>
</html>
	"""
	print >>HTML_OUT,html_subfix
	HTML_OUT.close()

	gj.printFuncRun('parse_html')

def rect_resize(coord_ls, plus=40):
	point1 = coord_ls[0:2]
	point2 = coord_ls[2:4]
	angle = math.atan(abs(point1[1]-point2[1])/float(abs(point1[0]-point2[0]))) if point1[0] != point2[0] else 0
	coord_ls_new = []
	coord_ls_new.append(point1[0]-plus*math.sin(angle))
	coord_ls_new.append(point1[1]+plus*math.cos(angle))
	coord_ls_new.append(point1[0]+plus*math.sin(angle))
	coord_ls_new.append(point1[1]-plus*math.cos(angle))
	coord_ls_new.append(point2[0]-plus*math.sin(angle))
	coord_ls_new.append(point2[1]+plus*math.cos(angle))
	coord_ls_new.append(point2[0]+plus*math.sin(angle))
	coord_ls_new.append(point2[1]-plus*math.cos(angle))
	return coord_ls_new


def coordinate_resize(coord_ls, y_plus=2, bias_range=2):
	print "coord_ls",coord_ls
	point_num = len(coord_ls)/2
	coord_ls_new = []
	for n in xrange(point_num):
		if n == 0:
			coord_ls_new.append(coord_ls[2*n])
			coord_ls_new.append(coord_ls[2*n+1])
			continue
		if coord_ls[2*n] == coord_ls_new[-2] and coord_ls[2*n+1] == coord_ls_new[-1]:
			continue
		else:
			coord_ls_new.append(coord_ls[2*n])
			coord_ls_new.append(coord_ls[2*n+1])
	print "coord_ls_new",coord_ls_new

	point_pair_num = point_num/2
	pair_success_num = 0
	for n in xrange(point_pair_num):
		point1 = coord_ls_new[n*2:n*2+2]
		point2 = coord_ls_new[::-1][n*2:n*2+2]
		if point1[0] == point2[1] and point1[1] == point2[0]:
			pair_success_num += 1
			print "[paired]",point1,point2
		elif abs(point1[0]-point2[1]) <= bias_range and abs(point1[1] - point2[0]) <= bias_range:
			pair_success_num += 1
			print "[paired]",point1,point2,"[within bias range]"
		else:
			print "[NOT paired]",point1,point2
	if pair_success_num/float(point_pair_num) > 0.8:
		print "[Line]",
		shape = "line"
	else:
		print "[polygon]",
		shape = 'polygon'

	if shape == 'polygon':
		print 
		return coord_ls
	else:
		coord_ls_resize = []
		for n in xrange(len(coord_ls_new)/2):
			if n <= point_pair_num:
				coord_ls_resize.append(coord_ls_new[2*n])
				coord_ls_resize.append(coord_ls_new[2*n+1]-y_plus)
			else:
				coord_ls_resize.append(coord_ls_new[2*n])
				coord_ls_resize.append(coord_ls_new[2*n+1]+y_plus)
		print coord_ls_resize
		return coord_ls_resize
	#sys.exit()
	return

def main():
    fetch_entry_data(mode='gene_id', id='ENSG00000100320')


if __name__ == "__main__":
    main()
