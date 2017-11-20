#!/bin/bash

### This is Anno_pipe, an easy and accessible pipeline
### For annotation of large datasets of proteins
### Author: Célio Dias Santos Júnior (celio.diasjunior@gmail.com)
### Usage: 
### Make an script qsub with SGE chunk and call
### 
###
###########################################################
### SGE Chunk for servers #################################
###########################################################
###$ -S /bin/bash
###$ -o /path/to/msg/REPORT_$JOB_ID #### change to your folder 
###$ -j y
###$ -cwd
###$ -V
###$ -pe qiime $3 ### change the number of processors
###$ -P protist ### change your group
###$ -N annotation
###$ -q all.q,share.q
##
## anno_pipe.sh <input_protein_fasta> <output_folder> <number_of_processors>
###########################################################

###########################################################
### Verifying #############################################
###########################################################

cpu=$3

if [ -z $1 ]; 
then 
	echo "## Usage:
## anno_pipe.sh <input_fasta_file> <output_folder> <threads>
[STATUS .. There is no input_file    ]
[STATUS .. input must be a fasta     ]
[STATUS .. review entering code //   ]"
	exit
elif [ -z $2 ]; 
then
	echo "## Usage:
## anno_pipe.sh <input_fasta_file> <output_folder> <threads>
[STATUS .. There is no output_folder]
[STATUS .. Output_folder don't exist]
[STATUS .. review entering code //  ]"
	exit
elif [ -z $3 ];
then
	echo "## Usage:
## anno_pipe.sh <input_fasta_file> <output_folder> <threads>
[STATUS .. No threads were associated]
[STATUS .. Assuming threads=4        ]
[STATUS .. Continuing                ]"
	cpu=4
fi
################################


###########################################################
### Setting up adresses ###################################
###########################################################

echo "H++ :: Setting up the libraries"

### For bioclusters users
### Decompress my previous made diamond databases
### /scratch2/cdias/databases/diamond_databases.tar.gz into
### your working folder and set these addresses

diamond=/share/apps/diamond/0.9.5/diamond
KEGG=/path/to/KEGG/diamond/database
COG=/path/to/COG/diamond/database
Uniprot=/path/to/Uniprot/diamond/database
CAMERA=/path/to/CAMERA/diamond/database
hmmer=/share/apps/hmmer/hmmer_latest/binaries/hmmsearch
EGGNOG=/gendata/databases/eggNOG/eggNOG_latest/NOG.hmm
PFAM=/gendata/databases/Pfam/Pfam_latest/Pfam-A.hmm

### Auxilliary scripts here 

hmm_out=/path_to/parse_hmmout.py
hmm_out2=/path_to/parse_egg_bdg.py
ref_folder=/path_to/ref_folder

###########################################################
### Starting processes ####################################
###########################################################

echo "H++ :: PFAM Searching"

$hmmer -o $2/PFAM_out --domtblout $2/PFAM_anno.domtblout --pfamtblout $2/PFAM_anno.pfamtblout -E 1e-3 --cpu $cpu $PFAM $1

echo "H++ :: PFAM Parsing"

python2 $hmm_out -i $2/PFAM_anno.domtblout -e 1e-5 --overlap 0 > $2/PFAM.parsed.1

cat $2/PFAM.parsed.1 | sort -k1,1 | uniq | awk -F'\t' -v OFS=',' '{x=$1;$1="";a[x]=a[x]$0}END{for(x in a)print x,a[x]}' | sed 's/,,/\t/g' > $2/PFAM.parsed.2

sort -k1,1 $2/PFAM.parsed.2 > $2/PFAM.parsed.25; rm -rf $2/PFAM.parsed.2; mv $2/PFAM.parsed.25 $2/PFAM.parsed.2

join $2/PFAM.parsed.2 $ref_folder/all.hmm.desc.sorted > $2/PFAM.functions; rm -rf PFAM.parsed*

echo "H++ :: EGGNOG Searching"

$hmmer -o $2/EggNOG_out --domtblout $2/EGGNOG_anno.domtblout --pfamtblout $2/EGGNOG_anno.pfamtblout -E 1e-3 --cpu $cpu $EGGNOG $1 

echo "H++ :: EGGNOG Parsing"

FILESIZE=$(stat -c%s "$2/EGGNOG_anno.domtblout")

if (( $FILESIZE < 2050000000 ))
then
	python2 $hmm_out2 -i $2/EGGNOG_anno.domtblout -e 1e-5 --overlap 0 > $2/EGG.parsed.1

	cat $2/EGG.parsed.1 | sort -k1,1 | uniq | awk -F'\t' -v OFS=',' '{x=$1;$1="";a[x]=a[x]$0}END{for(x in a)print x,a[x]}' | sed 's/,,/\t/g' > $2/EGG.parsed.2
else
	echo "H++ :: Your EGGNOG result was greater than 2GB..."

	echo "H++ :: To avoid problems we are skipping traditional parsing"

	echo "H++ :: Parsing in alternate way... Initiating Shell approach"

	cat $2/EGGNOG_anno.domtblout | grep -v '^#' | awk '{print $1,$3,$4,$6,$13,$16,$17,$18,$19}' | sed 's/ /\t/g' | \

	sort -k 3,3 -k 8n -k 9n | \

	perl -e 'while(<>){chomp;@a=split;next if $a[-1]==$a[-2];push(@{$b{$a[2]}},$_);}foreach(sort keys %b){@a=@{$b{$_}};for($i=0;$i<$#a;$i++){@b=split(/\t/,$a[$i]);@c=split(/\t/,$a[$i+1]);$len1=$b[-1]-$b[-2];$len2=$c[-1]-$c[-2];$len3=$b[-1]-$c[-2];if($len3>0 and ($len3/$len1>0.5 or $len3/$len2>0.5)){if($b[4]<$c[4]){splice(@a,$i+1,1);}else{splice(@a,$i,1);}$i=$i-1;}}foreach(@a){print $_."\n";}}' | \
	perl -e 'while(<>){chomp;@a=split(/\t/,$_);if(($a[-1]-$a[-2])>80){print $_,"\t",($a[-3]-$a[-4])/$a[1],"\n" if $a[4]<1e-5;}else{print $_,"\t",($a[-3]-$a[-4])/$a[1],"\n" if $a[4]<1e-3;}}' | awk '$NF>0.3' | sort -k 3 -k 8,9g > EGGNOG.parsed.final.tbl

fi

echo "H++ :: KEGG Searching"

$diamond blastp -d $KEGG --id 25.0 --query-cover 50.0 -q $1 -o $2/KEGG_diamond.tbl -f 6 -p $cpu

echo "H++ :: KEGG Parsing"

LANG=C; LC_ALL=C; sort -k1,1 -k12,12gr -k11,11g -k3,3gr $2/KEGG_diamond.tbl | sort -u -k1,1 --merge > $2/KEGG_filtered.tbl

cut -f 1-3 $2/KEGG_filtered.tbl > $2/starting

cut -f 2 $2/starting > $2/target

cut -f 1,3 $2/starting > $2/query

echo "H++ :: KEGG parsing... 1st stage completed"

paste -d '\t' $2/target $2/query > $2/file.1

sort -k1,1 $2/file.1 > $2/file.2

rm -rf $2/query $2/target $2/starting $2/file.1

echo "H++ :: KEGG parsing... 2nd stage completed"

join $2/file.2 $ref_folder/ko.sorted > $2/KEGG_1; rm -rvf $2/file.2

echo "H++ :: KEGG parsing... 3rd stage completed"

sed -i 's/ /\t/g' $2/KEGG_1

cut -f2,4 $2/KEGG_1 > $2/KEGG_2

cat $2/KEGG_2 | sort -k1,1 | uniq | awk -F'\t' -v OFS=',' '{x=$1;$1="";a[x]=a[x]$0}END{for(x in a)print x,a[x]}' | sed 's/,,/\t/g' > $2/KEGG_GO

echo "H++ :: KEGG parsing... stage completed"

echo "H++ :: CAMERA Searching"

$diamond blastp -d $CAMERA --id 25.0 --query-cover 50.0 -q $1 -o $2/CAMERA_diamond.tbl -f 6 -p $cpu

echo "H++ :: CAMERA Parsing"

LANG=C; LC_ALL=C; sort -k1,1 -k12,12gr -k11,11g -k3,3gr $2/CAMERA_diamond.tbl | sort -u -k1,1 --merge > $2/CAMERA_filtered.tbl

echo "H++ :: COG Searching"

$diamond blastp -d $COG --id 25.0 --query-cover 50.0 -q $1 -o $2/COG_diamond.tbl -f 6 -p $cpu

echo "H++ :: COG Parsing"

LANG=C; LC_ALL=C; sort -k1,1 -k12,12gr -k11,11g -k3,3gr $2/COG_diamond.tbl | sort -u -k1,1 --merge > $2/COG_filtered.tbl

cut -f 1-3 $2/COG_filtered.tbl > $2/starting

cut -f 2 $2/starting > $2/target

sed -i 's/gi|//g' $2/target

sed -i 's/|.*//g' $2/target

cut -f 1,3 $2/starting > $2/query

echo "H++ :: COG parsing... 1st stage completed"

paste -d '\t' $2/target $2/query > $2/file.1

sort -k1,1 $2/file.1 > $2/file.2

rm -rf $2/query $2/target $2/starting $2/file.1

echo "H++ :: COG parsing... 2nd stage completed"

join $2/file.2 $ref_folder/list_COG > $2/COG.1; rm -rf $2/file.2

echo "H++ :: COG parsing... 3rd stage completed"

sed -i 's/ /\t/g' $2/COG.1

cut -f2,4 $2/COG.1 > $2/COG.2

cat $2/COG.2 | sort -k1,1 | uniq | awk -F'\t' -v OFS=',' '{x=$1;$1="";a[x]=a[x]$0}END{for(x in a)print x,a[x]}' | sed 's/,,/\t/g' > $2/COG_results

echo "H++ :: COG parsing... stage completed"

echo "H++ :: Uniprot Searching"

$diamond blastp -d $Uniprot --id 25.0 --query-cover 50.0 -q $1 -o $2/UNIPROT_diamond.tbl -f 6 -p $cpu

echo "H++ :: Uniprot Parsing"

LANG=C; LC_ALL=C; sort -k1,1 -k12,12gr -k11,11g -k3,3gr $2/UNIPROT_diamond.tbl | sort -u -k1,1 --merge > $2/Uniprot_filtered.tbl

cut -f 1-3 $2/Uniprot_filtered.tbl > $2/starting

cut -f 2 $2/starting > $2/target

cut -f 1,3 $2/starting > $2/query

echo "H++ :: Uniprot parsing... 1st stage completed"

paste -d '\t' $2/target $2/query > $2/file.1

sort -k1,1 $2/file.1 > $2/file.2

rm -rf $2/query $2/target $2/starting $2/file.1

echo "H++ :: Uniprot parsing... 2nd stage completed"

join $2/file.2 $ref_folder/id_map_uniprot_cog > $2/Uniprot_1; rm -rvf $2/file.2

echo "H++ :: Uniprot parsing... 3rd stage completed"

sed -i 's/ /\t/g' $2/Uniprot_1

cut -f2,4 $2/Uniprot_1 > $2/Uniprot_2

cat $2/Uniprot_2 | sort -k1,1 | uniq | awk -F'\t' -v OFS=',' '{x=$1;$1="";a[x]=a[x]$0}END{for(x in a)print x,a[x]}' | sed 's/,,/\t/g' > $2/Uniprot_GO

cut -f2 $2/Uniprot_1 > $2/ctgs

cut -f1 $2/Uniprot_1 > $2/names

paste -d'\t' $2/ctgs $2/names > $2/u; rm -rf $2/ctgs $2/names Uniprot_1 Uniprot_2

sort -k1,1 $2/u | uniq | awk -F'\t' -v OFS=',' '{x=$1;$1="";a[x]=a[x]$0}END{for(x in a)print x,a[x]}' | sed 's/,,/\t/g' > $2/Uniprot_homologues; rm -rf $2/u

echo "H++ :: Uniprot parsing... stage completed"
