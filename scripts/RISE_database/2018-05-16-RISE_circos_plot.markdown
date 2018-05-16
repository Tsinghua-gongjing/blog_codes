---
layout: post
category: "genomics"
title:  "How to plot circos graph in RISE database"
tags: [genomics, RISE]
---

[Circos plot](http://circos.ca/) has been widely used to display genomic element interactions along with annotations. In [RISE database](http://rise.life.tsinghua.edu.cn/index.html) we also apply visualise RNA-RNA interactions in this format.

![img](http://rise.life.tsinghua.edu.cn/circos-data/human/ENSG00000100320.convert.svg)

Here is how we retrieve data from mysql database, process the data into specific format and plot with Circos tool.


####  login to mysql and check the database information:

```
[zhangqf5@loginview02 bin]$ ./mysql -h bio04 -uzhangqf5 
Warning: Using a password on the command line interface can be insecure.
Welcome to the MySQL monitor.  Commands end with ; or \g.
Your MySQL connection id is 1
Server version: 5.6.34-log MySQL Community Server (GPL)

Copyright (c) 2000, 2016, Oracle and/or its affiliates. All rights reserved.

Oracle is a registered trademark of Oracle Corporation and/or its
affiliates. Other names may be trademarks of their respective
owners.

Type 'help;' or '\h' for help. Type '\c' to clear the current input statement.

```

#### list available databases:

```
mysql> show DATABASES;
+--------------------+
| Database           |
+--------------------+
| information_schema |
| RISE               |
| test               |
+--------------------+
3 rows in set (0.24 sec)
```

#### use RISE database and list tables:

```
mysql> use RISE;
Reading table information for completion of table and column names
You can turn off this feature to get a quicker startup with -A

Database changed
mysql> show tables;
+-----------------------------------------+
| Tables_in_RISE                          |
+-----------------------------------------+
| auth_group                              |
| auth_group_permissions                  |
| auth_permission                         |
| auth_user                               |
| auth_user_groups                        |
| auth_user_user_permissions              |
| biomart_hg38                            |
| biomart_hg38_mm10_longest               |
| biomart_hg38_mm10_longest_chr           |
| biomart_mm10                            |
| conservation                            |
| datas_conservation                      |
| datas_rri                               |
| django_admin_log                        |
| django_content_type                     |
| django_migrations                       |
| django_session                          |
| expressionaltasCellLine_mouse           |
| gencode_biomart_hg38_mm10_yeast         |
| gencode_biomart_hg38_mm10_yeast_repeats |
| gencode_hg38_mm10                       |
| human_DARNED                            |
| human_RADAR                             |
| human_RMBaseNm                          |
| human_RMBaseOthers                      |
| human_RMBasePseudoU                     |
| human_RMBasem5C                         |
| human_RMBasem6A                         |
| human_clipdb                            |
| human_editing                           |
| human_modification                      |
| human_pancan                            |
| human_ucscSnp142                        |
| mouse_clipdb                            |
| mouse_editing                           |
| mouse_modification                      |
| mouse_ucscSnp142                        |
| proteinaltasCellLine_human              |
| rri                                     |
| rri_bk                                  |
+-----------------------------------------+
40 rows in set (0.00 sec)
```

#### columns of specific tables:

```
mysql> describe rri;
+-----------------------------+-------------+------+-----+---------+-------+
| Field                       | Type        | Null | Key | Default | Extra |
+-----------------------------+-------------+------+-----+---------+-------+
| index                       | bigint(20)  | YES  | MUL | NULL    |       |
| chr1                        | varchar(10) | YES  |     | NULL    |       |
| start1                      | varchar(10) | YES  |     | NULL    |       |
| end1                        | varchar(10) | YES  |     | NULL    |       |
| chr2                        | varchar(10) | YES  |     | NULL    |       |
| start2                      | varchar(10) | YES  |     | NULL    |       |
| end2                        | varchar(10) | YES  |     | NULL    |       |
| rise_id                     | varchar(13) | YES  |     | NULL    |       |
| score                       | varchar(15) | YES  |     | NULL    |       |
| strand1                     | varchar(1)  | YES  |     | NULL    |       |
| strand2                     | varchar(1)  | YES  |     | NULL    |       |
| ensembl_gene_id1            | varchar(40) | YES  |     | NULL    |       |
| ensembl_gene_name1          | varchar(40) | YES  |     | NULL    |       |
| ensembl_gene_id2            | varchar(40) | YES  |     | NULL    |       |
| ensembl_gene_name2          | varchar(40) | YES  |     | NULL    |       |
| ensembl_gene_type1          | varchar(50) | YES  |     | NULL    |       |
| ensembl_gene_type2          | varchar(50) | YES  |     | NULL    |       |
| ensembl_gene_abstract_type1 | varchar(20) | YES  |     | NULL    |       |
| ensembl_gene_abstract_type2 | varchar(20) | YES  |     | NULL    |       |
| method                      | varchar(15) | YES  |     | NULL    |       |
| species                     | varchar(20) | YES  |     | NULL    |       |
| cell_line                   | varchar(20) | YES  |     | NULL    |       |
| cell_line2                  | varchar(20) | YES  |     | NULL    |       |
| pubmed_id                   | text        | YES  |     | NULL    |       |
| reference                   | text        | YES  |     | NULL    |       |
| pos1                        | varchar(25) | YES  |     | NULL    |       |
| pos2                        | varchar(25) | YES  |     | NULL    |       |
| chr1_fill                   | varchar(10) | YES  |     | NULL    |       |
| start1_fill                 | varchar(10) | YES  |     | NULL    |       |
| end1_fill                   | varchar(10) | YES  |     | NULL    |       |
| chr2_fill                   | varchar(10) | YES  |     | NULL    |       |
| start2_fill                 | varchar(10) | YES  |     | NULL    |       |
| end2_fill                   | varchar(10) | YES  |     | NULL    |       |
| genomic_context1            | varchar(10) | YES  |     | NULL    |       |
| genomic_context2            | varchar(10) | YES  |     | NULL    |       |
| pos1_phastCons              | varchar(15) | YES  |     | NULL    |       |
| pos1_phyloP                 | varchar(15) | YES  |     | NULL    |       |
| pos2_phastCons              | varchar(15) | YES  |     | NULL    |       |
| pos2_phyloP                 | varchar(15) | YES  |     | NULL    |       |
+-----------------------------+-------------+------+-----+---------+-------+
39 rows in set (0.06 sec)
```

```
mysql> describe gencode_biomart_hg38_mm10_yeast_repeats;
+------------------+------------+------+-----+---------+-------+
| Field            | Type       | Null | Key | Default | Extra |
+------------------+------------+------+-----+---------+-------+
| index            | bigint(20) | YES  | MUL | NULL    |       |
| gene_id          | text       | YES  |     | NULL    |       |
| gene_name        | text       | YES  |     | NULL    |       |
| gene_start       | bigint(20) | YES  |     | NULL    |       |
| gene_end         | bigint(20) | YES  |     | NULL    |       |
| gene_type        | text       | YES  |     | NULL    |       |
| tx_id            | text       | YES  |     | NULL    |       |
| tx_start         | bigint(20) | YES  |     | NULL    |       |
| tx_end           | bigint(20) | YES  |     | NULL    |       |
| tx_type          | text       | YES  |     | NULL    |       |
| tx_abstract_type | text       | YES  |     | NULL    |       |
| sequence_type    | text       | YES  |     | NULL    |       |
| chr              | text       | YES  |     | NULL    |       |
| start            | bigint(20) | YES  |     | NULL    |       |
| end              | bigint(20) | YES  |     | NULL    |       |
+------------------+------------+------+-----+---------+-------+
15 rows in set (0.01 sec)
```

#### run script [RISE\_circos\_plot.py](https://github.com/Tsinghua-gongjing/blog_codes/tree/master/scripts/RISE_database) to generate the circos plot:

```
/Share/home/zhangqf/usr/anaconda/bin/python RISE_circos_plot.py
```

#### output files

```
[zhangqf5@loginview02 duplex_test]$ lt
total 11K
-rw-rw---- 1 zhangqf5 zhangqf 394K May 16 21:26 ENSG00000100320.convert.png
-rw-rw---- 1 zhangqf5 zhangqf 133K May 16 21:26 ENSG00000100320.convert.resize.html
-rw-rw---- 1 zhangqf5 zhangqf 135K May 16 21:26 ENSG00000100320.convert.svg
-rw-rw---- 1 zhangqf5 zhangqf 369K May 16 21:26 ENSG00000100320.log
-rw-rw---- 1 zhangqf5 zhangqf 258K May 16 21:26 ENSG00000100320.convert.html
-rw-rw---- 1 zhangqf5 zhangqf  42K May 16 21:26 ENSG00000100320.anno.convert
-rw-rw---- 1 zhangqf5 zhangqf  37K May 16 21:26 ENSG00000100320.genestructure.convert
-rw-rw---- 1 zhangqf5 zhangqf  17K May 16 21:26 ENSG00000100320.link.convert
-rw-rw---- 1 zhangqf5 zhangqf 5.3K May 16 21:26 ENSG00000100320.linkgenes.convert
-rw-rw---- 1 zhangqf5 zhangqf 2.8K May 16 21:26 ENSG00000100320.linkgenes.txt.convert
-rw-rw---- 1 zhangqf5 zhangqf 349K May 16 21:26 ENSG00000100320.png
-rw-rw---- 1 zhangqf5 zhangqf 135K May 16 21:26 ENSG00000100320.svg
-rw-rw---- 1 zhangqf5 zhangqf 232K May 16 21:26 ENSG00000100320.html
-rw-rw---- 1 zhangqf5 zhangqf  42K May 16 21:26 ENSG00000100320.anno
-rw-rw---- 1 zhangqf5 zhangqf  37K May 16 21:26 ENSG00000100320.genestructure
-rw-rw---- 1 zhangqf5 zhangqf  22K May 16 21:26 ENSG00000100320.genestructure.introns
-rw-rw---- 1 zhangqf5 zhangqf  19K May 16 21:26 ENSG00000100320.genestructure.introns.txt
-rw-rw---- 1 zhangqf5 zhangqf  17K May 16 21:26 ENSG00000100320.link
-rw-rw---- 1 zhangqf5 zhangqf 5.3K May 16 21:26 ENSG00000100320.linkgenes
-rw-rw---- 1 zhangqf5 zhangqf 2.8K May 16 21:26 ENSG00000100320.linkgenes.txt
-rw-rw---- 1 zhangqf5 zhangqf 1.9K May 16 21:26 ENSG00000100320.linkgenes.zooms
-rw-rw---- 1 zhangqf5 zhangqf 2.1K May 16 21:26 ENSG00000100320.linkgenes.zooms2
```

Note:

ENSG00000100320.convert.png/svg is the final plot correspond to the figure above or in the database for RNA [RBFOX2](http://rise.life.tsinghua.edu.cn/search-result.html?species=human&genename=rbfox2)

#### track files (before log2 transformed): desired file formats are [here ENSG00000100320(RBFOX2) as an example](https://github.com/Tsinghua-gongjing/blog_codes/tree/master/scripts/RISE_database/ENSG00000100320)

file gene structure:

```
[zhangqf5@loginview02 duplex_test]$ head -5 ENSG00000100320.genestructure
hs9 94175952 94176043 gene_id=ENSG00000199165,gene_name=MIRLET7A1,gene_type=miRNA,tx_id=ENST00000362295,sequence_type=Exon,id=ENSG00000199165|MIRLET7A1|94175952|94176043|miRNA|ENST00000362295|94175952|94176043|miRNA|miRNA|Exon|9|94175952|94176043
hs10 133246478 133247050 gene_id=ENSG00000166917,gene_name=MIR202HG,gene_type=ncRNA,tx_id=ENST00000553459,sequence_type=Exon,id=ENSG00000166917|MIR202HG|133246478|133247891|antisense|ENST00000553459|133246478|133247891|antisense|ncRNA|Exon|10|133246478|133247050
hs10 133247791 133247891 gene_id=ENSG00000166917,gene_name=MIR202HG,gene_type=ncRNA,tx_id=ENST00000553459,sequence_type=Exon,id=ENSG00000166917|MIR202HG|133246478|133247891|antisense|ENST00000553459|133246478|133247891|antisense|ncRNA|Exon|10|133247791|133247891
hs17 77319515 77320326 gene_id=ENSG00000184640,gene_name=SEPT9,gene_type=protein_coding,tx_id=ENST00000329047,sequence_type=5UTR,id=ENSG00000184640|SEPT9|77280569|77500596|protein_coding|ENST00000329047|77319515|77500592|protein_coding|protein_coding|5UTR|17|77319515|77320326
hs17 77498656 77500592 gene_id=ENSG00000184640,gene_name=SEPT9,gene_type=protein_coding,tx_id=ENST00000329047,sequence_type=3UTR,id=ENSG00000184640|SEPT9|77280569|77500596|protein_coding|ENST00000329047|77319515|77500592|protein_coding|protein_coding|3UTR|17|77498656|77500592
``` 

file annotations:

```
[zhangqf5@loginview02 duplex_test]$ head -5 ENSG00000100320.anno
hs22 35744177 35744178 1 anno_category=modification,anno_source=RMBase,anno_type=m6A,anno_gene_info=chr22|35744164|35744204|RISE0255307_2|2|-|ENSG00000100320|RBFOX2|protein_coding|protein_coding,anno_full=chr22|35744177|35744178|m6A_site_86471|0|-|m6A_site_86471|m6A|1|GSE60213|25799998|RBFOX2&RBFOX2|retained_intron&protein_coding|exon&utr3|GTGAGACCCCTGCAAATGGGACAGCCCCCCAGTTCATGAGG|modification|RMBase|https://www.ncbi.nlm.nih.gov/pubmed/25799998|ENSG00000100320|ENST00000405409|3UTR|1,id=chr22|35744164|35744204|RISE0255307_2|2|-|ENSG00000100320|RBFOX2|protein_coding|protein_coding|chr22|35744177|35744178|m6A_site_86471|0|-|m6A_site_86471|m6A|1|GSE60213|25799998|RBFOX2&RBFOX2|retained_intron&protein_coding|exon&utr3|GTGAGACCCCTGCAAATGGGACAGCCCCCCAGTTCATGAGG|modification|RMBase|https://www.ncbi.nlm.nih.gov/pubmed/25799998|ENSG00000100320|ENST00000405409|3UTR|1
hs22 35744192 35744193 1 anno_category=modification,anno_source=RMBase,anno_type=m6A,anno_gene_info=chr22|35744164|35744204|RISE0255307_2|2|-|ENSG00000100320|RBFOX2|protein_coding|protein_coding,anno_full=chr22|35744192|35744193|m6A_site_86472|0|-|m6A_site_86472|m6A|1|GSE60213|25799998|RBFOX2&RBFOX2|retained_intron&protein_coding|exon&utr3|CCCTACTGAAGTGACGTGAGACCCCTGCAAATGGGACAGCC|modification|RMBase|https://www.ncbi.nlm.nih.gov/pubmed/25799998|ENSG00000100320|ENST00000405409|3UTR|1,id=chr22|35744164|35744204|RISE0255307_2|2|-|ENSG00000100320|RBFOX2|protein_coding|protein_coding|chr22|35744192|35744193|m6A_site_86472|0|-|m6A_site_86472|m6A|1|GSE60213|25799998|RBFOX2&RBFOX2|retained_intron&protein_coding|exon&utr3|CCCTACTGAAGTGACGTGAGACCCCTGCAAATGGGACAGCC|modification|RMBase|https://www.ncbi.nlm.nih.gov/pubmed/25799998|ENSG00000100320|ENST00000405409|3UTR|1
hs5 6668743 6668744 1 anno_category=RBP,anno_source=CLIPDB,anno_type=RBP,anno_gene_info=chr5|6668720|6668754|RISE0224410_1|2|+|ENSG00000145545|SRD5A1|protein_coding|protein_coding,anno_full=chr5|6668743|6668744|human_RBP_CLIPdb_iCLIP_CIMS_2234997|0|+|PTBP1|iCLIP&CIMS|HEK293T|GSE57278&GSM1378377|9|25219497|https://www.ncbi.nlm.nih.gov/pubmed/25219497|ENSG00000145545|ENST00000510531|3UTR|1|iCLIP|CIMS|human,id=chr5|6668720|6668754|RISE0224410_1|2|+|ENSG00000145545|SRD5A1|protein_coding|protein_coding|chr5|6668743|6668744|human_RBP_CLIPdb_iCLIP_CIMS_2234997|0|+|PTBP1|iCLIP&CIMS|HEK293T|GSE57278&GSM1378377|9|25219497|https://www.ncbi.nlm.nih.gov/pubmed/25219497|ENSG00000145545|ENST00000510531|3UTR|1|iCLIP|CIMS|human
hs5 6668704 6668725 1 anno_category=RBP,anno_source=CLIPDB,anno_type=RBP,anno_gene_info=chr5|6668720|6668754|RISE0224410_1|2|+|ENSG00000145545|SRD5A1|protein_coding|protein_coding,anno_full=chr5|6668704|6668725|human_RBP_CLIPdb_iCLIP_CIMS_4297059|0|+|TARDBP|iCLIP&CIMS|Brain|E-MTAB-530&ERR039849|4|21358640|https://www.ncbi.nlm.nih.gov/pubmed/21358640|ENSG00000145545|ENST00000510531|3UTR|5|iCLIP|CIMS|human,id=chr5|6668720|6668754|RISE0224410_1|2|+|ENSG00000145545|SRD5A1|protein_coding|protein_coding|chr5|6668704|6668725|human_RBP_CLIPdb_iCLIP_CIMS_4297059|0|+|TARDBP|iCLIP&CIMS|Brain|E-MTAB-530&ERR039849|4|21358640|https://www.ncbi.nlm.nih.gov/pubmed/21358640|ENSG00000145545|ENST00000510531|3UTR|5|iCLIP|CIMS|human
hs5 6668749 6668770 1 anno_category=RBP,anno_source=CLIPDB,anno_type=RBP,anno_gene_info=chr5|6668720|6668754|RISE0224410_1|2|+|ENSG00000145545|SRD5A1|protein_coding|protein_coding,anno_full=chr5|6668749|6668770|human_RBP_CLIPdb_iCLIP_CIMS_4435494|0|+|TIAL1|iCLIP&CIMS|HeLa|E-MTAB-432&ERR039780-ERR039783-ERR039784-ERR039785|7|21048981|https://www.ncbi.nlm.nih.gov/pubmed/21048981|ENSG00000145545|ENST00000510531|3UTR|5|iCLIP|CIMS|human,id=chr5|6668720|6668754|RISE0224410_1|2|+|ENSG00000145545|SRD5A1|protein_coding|protein_coding|chr5|6668749|6668770|human_RBP_CLIPdb_iCLIP_CIMS_4435494|0|+|TIAL1|iCLIP&CIMS|HeLa|E-MTAB-432&ERR039780-ERR039783-ERR039784-ERR039785|7|21048981|https://www.ncbi.nlm.nih.gov/pubmed/21048981|ENSG00000145545|ENST00000510531|3UTR|5|iCLIP|CIMS|human
```

file links (RRI):

```
[zhangqf5@loginview02 duplex_test]$ head -5 ENSG00000100320.link
hs5 6668720 6668754 hs22 35740046 35740082 category=GlobalAnalysis,method=PARIS,duplex=derived,gene_id1=ENSG00000145545,gene_name1=SRD5A1,gene_type1=protein_coding,gene_id2=ENSG00000100320,gene_name2=RBFOX2,gene_type2=protein_coding,id=chr5|6668720|6668754|chr22|35740046|35740082|RISE0224410|2|+|-|ENSG00000145545|SRD5A1|ENSG00000100320|RBFOX2|protein_coding|protein_coding|protein_coding|protein_coding|PARIS|human|Hela(highRNase)|Hela|27180905|https://www.ncbi.nlm.nih.gov/pubmed/27180905|chr5:6668720-6668754|chr22:35740046-35740082|chr5|6668720|6668754|chr22|35740046|35740082|3UTR|3UTR|0.000441176|-0.350088|0.0715833|0.0108889
hs22 35742350 35742375 hs17 41866974 41867001 category=GlobalAnalysis,method=PARIS,duplex=derived,gene_id1=ENSG00000100320,gene_name1=RBFOX2,gene_type1=protein_coding,gene_id2=ENSG00000131473,gene_name2=ACLY,gene_type2=protein_coding,id=chr22|35742350|35742375|chr17|41866974|41867001|RISE0233924|2|-|-|ENSG00000100320|RBFOX2|ENSG00000131473|ACLY|protein_coding|protein_coding|protein_coding|protein_coding|PARIS|human|Hela(highRNase)|Hela|27180905|https://www.ncbi.nlm.nih.gov/pubmed/27180905|chr22:35742350-35742375|chr17:41866974-41867001|chr22|35742350|35742375|chr17|41866974|41867001|3UTR|3UTR|0.2758|0.29404|0.733481|0.895222
hs22 35741848 35741867 hs8 101129130 101129157 category=GlobalAnalysis,method=PARIS,duplex=derived,gene_id1=ENSG00000100320,gene_name1=RBFOX2,gene_type1=protein_coding,gene_id2=ENSG00000253355,gene_name2=KB-1460A1,gene_type2=ncRNA,id=chr22|35741848|35741867|chr8|101129130|101129157|RISE0233925|2|-|+|ENSG00000100320|RBFOX2|ENSG00000253355|KB-1460A1|protein_coding|processed_transcript|protein_coding|ncRNA|PARIS|human|Hela(highRNase)|Hela|27180905|https://www.ncbi.nlm.nih.gov/pubmed/27180905|chr22:35741848-35741867|chr8:101129130-101129157|chr22|35741848|35741867|chr8|101129130|101129157|3UTR|Exon|0.000315789|-0.186737|0.0836296|0.0691481
hs22 35741792 35741815 hs6 123471015 123471034 category=GlobalAnalysis,method=PARIS,duplex=derived,gene_id1=ENSG00000100320,gene_name1=RBFOX2,gene_type1=protein_coding,gene_id2=ENSG00000235535,gene_name2=RP11-532N4,gene_type2=ncRNA,id=chr22|35741792|35741815|chr6|123471015|123471034|RISE0233926|2|-|+|ENSG00000100320|RBFOX2|ENSG00000235535|RP11-532N4|protein_coding|antisense|protein_coding|ncRNA|PARIS|human|Hela(highRNase)|Hela|27180905|https://www.ncbi.nlm.nih.gov/pubmed/27180905|chr22:35741792-35741815|chr6:123471015-123471034|chr22|35741792|35741815|chr6|123471015|123471034|3UTR|Exon|0.00156522|0.232348|0.0652632|0.296211
hs22 35740318 35740336 hs2 97735297 97735318 category=GlobalAnalysis,method=PARIS,duplex=derived,gene_id1=ENSG00000100320,gene_name1=RBFOX2,gene_type1=protein_coding,gene_id2=ENSG00000115085,gene_name2=ZAP70,gene_type2=protein_coding,id=chr22|35740318|35740336|chr2|97735297|97735318|RISE0233927|2|-|+|ENSG00000100320|RBFOX2|ENSG00000115085|ZAP70|protein_coding|protein_coding|protein_coding|protein_coding|PARIS|human|Hela(highRNase)|Hela|27180905|https://www.ncbi.nlm.nih.gov/pubmed/27180905|chr22:35740318-35740336|chr2:97735297-97735318|chr22|35740318|35740336|chr2|97735297|97735318|3UTR|CDS|0.895778|1.97189|0.903|4.62224
```

file link genes:

```
[zhangqf5@loginview02 duplex_test]$ head -5 ENSG00000100320.linkgenes
hs9 94175952 94176043 gene_id=ENSG00000199165,gene_name=MIRLET7A1,gene_type=miRNA,tx_id=ENST00000362295,id=ENSG00000199165|MIRLET7A1|94175952|94176043|miRNA|ENST00000362295|94175952|94176043|miRNA|miRNA|Exon|9|94175952|94176043
hs17 41866908 41918987 gene_id=ENSG00000131473,gene_name=ACLY,gene_type=protein_coding,tx_id=ENST00000352035,id=ENSG00000131473|ACLY|41866908|41930542|protein_coding|ENST00000352035|41866908|41918987|protein_coding|protein_coding|5UTR|17|41913874|41913896
hs10 133246478 133247891 gene_id=ENSG00000166917,gene_name=MIR202HG,gene_type=miRNA,tx_id=ENST00000553459,id=ENSG00000166917|MIR202HG|133246478|133247891|antisense|ENST00000553459|133246478|133247891|antisense|ncRNA|Exon|10|133246478|133247050
hs1 629640 630683 gene_id=ENSG00000225630,gene_name=MTND2P28,gene_type=pseudogene,tx_id=ENST00000457540,id=ENSG00000225630|MTND2P28|629640|630683|unprocessed_pseudogene|ENST00000457540|629640|630683|unprocessed_pseudogene|pseudogene|Exon|1|629640|630683
hsX 53556223 53556341 gene_id=ENSG00000271886,gene_name=MIR98,gene_type=miRNA,tx_id=ENST00000606724,id=ENSG00000271886|MIR98|53556223|53556341|miRNA|ENST00000606724|53556223|53556341|miRNA|miRNA|Exon|X|53556223|53556341
```

#### generated .conf file, which control the layout of whole plot, e.g., [ENSG00000100320.conf](https://github.com/Tsinghua-gongjing/blog_codes/blob/master/scripts/RISE_database/ENSG00000100320/ENSG00000100320.conf)
