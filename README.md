# Document-description
tool usage and experimental data for master's work
## 1.Tool usage
### 1.1 self error correction tool
**LCAT**  
function: long-read error correction algorithm for transcriptome sequencing data  
website: https://github.com/luckylyw/LCAT  
usage: see at github-README.md  
**MECAT**<br>
function: an ultra-fast Mapping, Error Correction and de novo Assembly Tools for single molecula sequencing (SMRT) reads.<br>
website: https://github.com/xiaochuanle/MECAT<br>
usage:<br>
```
mecat2pw -j 0 -d SRR1284073.fastq -o SRR1284073.fastq.pm.can -w wrk_dir -t 16 7
mecat2cns -i 0 -t 16 SRR1284073.fastq.pm.can SRR1284073.fastq corrected_reads.fastq
```
**HALS**<br>
function:HALS is software that makes self-correction for long reads with fast speed and high throughput.<br>
website:https://github.com/xief001/hals<br>
install:<br>
```
git clone https://github.com/xief001/HALS.git
cd HALS/
make
vi ~/.bashrc
export PATH=/home/luoluo/tool/HALS/bin:$PATH
(修改runHALS.py中的路径)
```
usage:<br>
```
awk '{if(NR%4==1){print ">" substr($0, 2)}}{if(NR%4 == 2){print}}' ERR968960.fastq > ERR968960.fasta
python3 ~/tool/HALS/runHALS.py ERR1655118.fasta 1>hals.log 2>hals.err
```
**Canu**<br>
function:Canu is a fork of the Celera Assembler, designed for high-noise single-molecule sequencing<br>
website:https://github.com/marbl/canu<br>
install:<br>
```
git clone https://github.com/marbl/canu.git
cd canu/src
make -j <number of threads>
设置环境变量
安装 gnuplot-5.0.5
$ tar -zxvf gnuplot-5.0.5.tar.gz
$ configure --prefix=/public1/home/Serenity/installed_software/gnuplot-5.0.5
$ make
$ make install
```
usage:<br>
```
canu -correct  -p ecoli1 -d ecoli1  genomeSize=4.5m  -pacbio-raw  xx.fastq
canu -trim -p ecoli -d ecoli genomeSIze=4.8m -pacbio-corrected ecoli/ecoli.correctedReads.fasta.gz
canu -correct -p ecoli -d ecoli genomeSIze=4.5m -nanopore_raw oxford.fasta
canu -correct -p output.canu -d $workdir genomeSize=190136026 useGrid=false ovsMethod=sequential corOutCoverage=all minReadLength=100 minOverlapLength=100 -maxThreads=32 -maxMemory=120g corOverlapper=minimap -nanopore-raw $NANOPORE_READS
```
**FALCON**<br>
function:a set of tools for fast aligning long reads for consensus and assembly<br>
website:https://github.com/PacificBiosciences/FALCON/<br>
install:<br>
```curl -O https://downloads.pacbcloud.com/public/falcon/falcon-2018.08.08-21.41-py2.7-ucs4-beta.tar.gz
export LD_LIBRARY_PATH=/home/luoluo/tool/FALCON/lib:${LD_LIBRARY_PATH}
取消注释export PYTHONPATH=$PYTHONPATH:/home/luoluo/tool/FALCON/lib/python2.7/site-packages
export PATH=/home/luoluo/.nimble/bin:$PATH
export PATH=/home/luoluo/tool/FALCON/bin:$PATH
```
usage:<br>
```
awk '{if(NR%4==1){print ">" substr($0, 2)}}{if(NR%4 == 2){print}}' ERR968960.fastq > ERR968960.fasta
perl falcon_name_fasta.pl  --infile <fasta input>  --outfile  <fasta output>
复制并修改input.fofn
cp /home/luoluo/data/flacon/SRR1284073/fc_run.cfg .
看情况修改fc_run.cfg文件
fc_run.py fc_run.cfg 1>SRR7533629.1.log 2>SRR7533629.1.err &
tail -f SRR7533629.1.err
```
**LoRMA**<br>
function:accurate self-correction of errors in long reads using de Bruijn graphs<br>
website:https://www.cs.helsinki.fi/u/lmsalmel/LoRMA/<br>
usage:<br>
```
conda install -c bioconda lorma
lorma.sh ana.fasta -threads 20
```
### 1.2 reference-based long read error correction tool for RNA
**TranscriptClean**<br>
function:corrects mismatches, microindels, and noncanonical splice junctions in long reads that have been mapped to the genome. It is designed for use with sam files from the PacBio Iso-seq and Oxford Nanopore transcriptome sequencing technologies.<br>
website:https://github.com/mortazavilab/TranscriptClean<br>
usage:<br>
```cd ~/tool/TranscriptClean-master/
luoluo@bio10:~/tool/TranscriptClean-master$ python2 TranscriptClean.py --sam ~/TEST_DATA/.../pse_aln.sam --genome ~/TEST_DATA/../pse_ref.fasta --outprefix ~/TEST_DATA/.../result
```
**TLCR**<br>
function:transcriptome long-read error correction based on reference<br>
install:<br>
```
install TranscriptClean
```
usage:<br>
```
step 1: alignment raw reads to genome reference
minimap2 -ax splice -uf -C5 --MD ref.fa raw_reads.fastq > aln.sam
step 2:output pse_ref.fasta
python pse_ref.py aln.sam ref.fa raw_reads.fastq
step 3:correct error using TranscriptClean
cd ~/tool/TranscriptClean-master/
python2 TranscriptClean.py --sam ~/TEST_DATA/.../pse_aln.sam --genome ~/TEST_DATA/../pse_ref.fasta --outprefix ~/TEST_DATA/.../result
step 4:union result
python find_unmap_fa_sam.py pse_aln.sam SRR1163655.fasta raw_aln.sam  1 > find.log
```
### 1.3 align tool
**Minimap2**<br>
function:a versatile sequence alignment program that aligns DNA or mRNA sequences against a large reference database<br>
website:https://github.com/lh3/minimap2<br>
install:<br>
```
git clone https://github.com/lh3/minimap2
cd minimap2 && make
export PATH=/home/luoluo/tool/minimap2:$PATH
```
usage:<br>
```
DNA:
minimap2 -d NC_000913.3.fasta corrected_SRR7533629.fasta > corrected_to_ref.sam
RNA:
minimap2 -ax splice --eqx --MD -uf -C5 Mus_musculus.GRCm38.dna.primary_assembly.fa ERR2401483_proccessed_normalid.fasta > raw_md_aln.sam
```
### 1.4 evaluation tool
**FASTQC**<br>
```
fastqc -o wrk -t 20 SRR1284073.fastq
```
**samtools**<br>
```
samtools view -bS aln.sam > aln.ban
samtools sort aln.bam > aln_sorted.bam
samtools index aln_sorted.bam
```
**alignQC**<br>
function: evaluation<br>
website:[https://github.com/jason-weirather/AlignQC/wiki](https://github.com/jason-weirather/AlignQC<br>
install:<br>
```
$ git clone https://github.com/jason-weirather/AlignQC.git
$ cd AlignQC
$ pip install AlignQC
```
usage:<br>
```
alignqc analyze long_reads.bam -g ref_genome.fa -t ref_transcriptome.gtf -a genepred -o alignqc_report.xhtml
--no_genome --no_transcriptome //--no_annotation
```
**LR_EC_analyser**<br>
function:LR_EC_analyser stands for Long Read Error Correction analyser.<br>
website:https://github.com/leoisl/LR_EC_analyser<br>
usage:<br>
```
cd ~/tool/LR_EC_analyser/
mkdir mouse
cp raw_aln.bam,mecat_aln.bam,my_aln_bam ./mouse
source venv/bin/activate
luoluo@bio10:~/tool/LR_EC_analyser/venv/bin$ls -lh
-rwxrwxr-x 1 luoluo luoluo   47 1月  10 16:23 orca
lrwxrwxrwx 1 luoluo luoluo   36 1月  10 16:47 orca_core ->
head orca.js
#!/bin/bash
    xvfb-run -a -s '-screen 0 640x480x24' venv/bin/orca_core "$@“
(venv) luoluo@bio10:~/tool/LR_EC_analyser$ python run_LR_EC_analyser.py --genome /home/luoluo/TEST_DATA/mouse/Mus_musculus.GRCm38.dna.primary_assembly.fa --gtf /home/luoluo/TEST_DATA/mouse/Mus_musculus.GRCm38.87.gtf  -t 20 -o mouse_pacbio/output --raw mouse_pacbio/raw_aln.bam --self mouse_pacbio/para2_aln.bam mouse_pacbio/my_aln.bam
python run_LR_EC_analyser.py --genome /home/luoluo/TEST_DATA/ana/GCF_000699085.1_ASM69908v1_genomic_anna.fasta --gtf /home/luoluo/TEST_DATA/ana/ref_gff/GCF_000699085.1_ASM69908v1_genomic.gtf  -t 20 -o ./anna_100/ --raw ./anna_100/raw_aln.bam --self ./anna_100/mecat_aln.bam ./anna_100/my_aln.bam
```
**seqkit**<br>
function:Sequential information extraction tool<br>
website:https://www.jianshu.com/p/471283080bd6<br>
install:<br>
```
conda install seqkit
```
usage:<br>
```
seqkit stats corrected_reads.fastq 1>seqkit_stats.log
file                               format  type  num_seqs     sum_len  min_len  avg_len  max_len
corrected_reads.fastq  FASTA   DNA     32,298  91,158,476    2,000  2,822.4    9,948
```
## 2.Experimental data
**ssh://luoluo:luoluo123@202.197.66.217:5009**<br>
### 2.1 LCAT-self error correction tool
four RNA long read datasets (/home/luoluo/LCAT_data/)	<br>
* mouse_nanopore <br>
* human_nanopore <br>
* zebra_pacbio <br>
* anna_pacbio<br>
In the folder corresponding to each species data, there are the following contents: <br>
* genome_reference.fasta <br>
* annotation.gtf.gz <br>
* raw_read.fasta <br>
* xx_lr_ec result <br>
* dir"raw/mecat/lcat"<br>
In the folder "raw/mecat/lcat", there are the following contents:<br>

content | file
---- | ----- 
alignment file | xxx.pm.can/candidatex.txt
correction result | corrected_reads.fastq
map result to ref | aln*
result read information | seqkit_stats.log
### 2.2 TLCR reference-guided error correction tool
three RNA long read datasets (/home/luoluo/TLCR_data):<br>
* zebra_pacbio<br>
* anna_pacbio<br>
* human_pacbio<br>
In the folder corresponding to each species data(zebra/human/anna), there are the following contents:<br>

content | file
---- | ----- 
reference | genome_ref.fasta
annotation | annotation.gtf.gz
raw reads | raw_read.fasta
tc result | rc.fa
tlcr result | union_result.fa
