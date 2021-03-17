#!/bin/bash

## 	Author: Laurens Lambrechts
##	Contact: laurens.lambrechts@ugent.be
## 	Date: September 2020

#set working directory
wrk=~/Desktop/HIV_STAR_assembly

#set RunID (input files should be placed in folder with this name)
run=STAR

#define project folder
prj=$wrk/$run

#set database folder
db=/$wrk/db/

#set virus for filtering
virus=HIV
blast="$virus"_refs

#set compute resources where m=GB t=CPU
m=30
t=10

############# Trimming #############

#change into project folder
cd $prj

#create folder to hold reads for processing (original will remain as backup)
mkdir $prj/1_reads


#quality trim reads with bbduk
while read i; do
bbduk.sh -Xmx"$m"g threads="$t" in=$prj/1_reads/"$i"_R1.fastq.gz in2=$prj/1_reads/"$i"_R2.fastq.gz out=stdout.fq entropy=0.7 | bbduk.sh threads="$t" interleaved=true in=stdin.fq out=stdout.fq ref=phix k=31 hdist=1 | bbduk.sh threads="$t" interleaved=true in=stdin.fq out=stdout.fq ref=adapters ktrim=r k=23 mink=7 | bbduk.sh threads="$t" interleaved=true in=stdin.fq out=$prj/1_reads/"$i"_trim_1.fastq.gz out2=$prj/1_reads/"$i"_trim_2.fastq.gz forcetrimleft=14 forcetrimright2=1 minlen=75 qtrim=rl trimq=25
done < $prj/IDs.list

#normalise read coverage for denovo assembly
while read i; do
bbnorm.sh -Xmx"$m"g threads="$t" in=$prj/1_reads/"$i"_trim_1.fastq.gz in2=$prj/1_reads/"$i"_trim_2.fastq.gz out=$prj/1_reads/"$i"_trim_norm_1.fastq.gz out2=$prj/1_reads/"$i"_trim_norm_2.fastq.gz target=100;
done < $prj/IDs.list 


#create folder for 2_ref_maps
mkdir $prj/2_ref_map

#map trim reads to reference genomes
while read i; do
bbmap.sh -Xmx"$m"g threads="$t" in=$prj/1_reads/"$i"_trim_1.fastq.gz in2=$prj/1_reads/"$i"_trim_2.fastq.gz outm=$prj/2_ref_map/"$i".ref_mapped.bam ref="$db"/blast/HXB2.fasta;
samtools sort -@ "$t" -o $prj/2_ref_map/"$i".ref_mapped.sorted.bam $prj/2_ref_map/"$i".ref_mapped.bam
rm $prj/2_ref_map/"$i".ref_mapped.bam
done < $prj/IDs.list 



#generate coverage maps
cd $prj/2_ref_map
while read i; do
qualimap bamqc -bam "$i".ref_mapped.sorted.bam  -outformat PDF -outfile "$i".ref_mapped.coverage.pdf --java-mem-size="$m"G
mv $prj/2_ref_map/"$i".ref_mapped.sorted_stats/"$i".ref_mapped.coverage.pdf $prj/2_ref_map/
rm -r $prj/2_ref_map/"$i".ref_mapped.sorted_stats
done < $prj/IDs.list
cd $prj

############# DE NOVO ASSEMBLY #############


#create folder for de novo assembly
mkdir $prj/3_contigs

#denovo assemble trimmed reads with megahit/1.2.9 and rename output contigs with library name
while read i; do
megahit --k-min 45 --k-max 95 --k-step 10 -t "$t" -1 $prj/1_reads/"$i"_trim_norm_1.fastq.gz -2 $prj/1_reads/"$i"_trim_norm_2.fastq.gz -o $prj/3_contigs/megahit_"$i" --out-prefix "$i" --min-contig-len 500;
sed "s/>k/>Draft_"$i"_k/g" $prj/3_contigs/megahit_"$i"/"$i".contigs.fa > $prj/3_contigs/megahit_"$i"/"$i".contigs.megahit.fa;
done < $prj/IDs.list

#export denovo assembly contig stats
while read i; do
grep ">" $prj/3_contigs/megahit_"$i"/"$i".contigs.megahit.fa | cut -d ">" -f2 | tr ' ' '\t' | sed 's/flag\=//g' | sed 's/multi\=//g' | sed 's/len\=//g' > $prj/3_contigs/megahit_"$i"/"$i".contigs.megahit.stats
done < $prj/IDs.list

#extract high coverage denovo assembly contigs
while read i; do
awk -F$'\t' '{OFS=FS}{if ($3>50) print $1}' $prj/3_contigs/megahit_"$i"/"$i".contigs.megahit.stats > $prj/3_contigs/megahit_"$i"/"$i".contigs.megahit.hicov.list
seqtk subseq $prj/3_contigs/megahit_"$i"/"$i".contigs.megahit.fa $prj/3_contigs/megahit_"$i"/"$i".contigs.megahit.hicov.list > $prj/3_contigs/megahit_"$i"/"$i".contigs.megahit.hicov.fa
done < $prj/IDs.list

#create folder for blasting
mkdir $prj/4_filter

#combine denovo contigs into single file
cat $prj/3_contigs/megahit_*/*.contigs.megahit.hicov.fa > $prj/4_filter/all_denovo.contigs.megahit.hicov.fasta


#blast to local HIV-1 database
export BLASTDB=$db/blast
blastn -query $prj/4_filter/all_denovo.contigs.megahit.hicov.fasta -db "$db"/blast/HXB2.fasta -evalue 1E-10 -num_threads "$t" -out $prj/4_filter/all_denovo.contigs.megahit.hicov.blastn_"$virus".txt -outfmt "6 qseqid qlen stitle sstart send pident length evalue sstrand"

#get top virus blast results for each contig
awk -F$'\t' '!seen[$1]++' $prj/4_filter/all_denovo.contigs.megahit.hicov.blastn_"$virus".txt > $prj/4_filter/all_denovo.contigs.megahit.hicov.blastn_"$virus".top.txt

#take column containing contig names
cut -f1 $prj/4_filter/all_denovo.contigs.megahit.hicov.blastn_"$virus".top.txt > $prj/4_filter/"$virus"_draft_genomes."$run".list 

#retrieve sequences
seqtk subseq $prj/4_filter/all_denovo.contigs.megahit.hicov.fasta $prj/4_filter/"$virus"_draft_genomes."$run".list > $prj/4_filter/"$virus"_draft_genomes."$run".fa

---


#align to reference
mafft --thread $t --reorder --adjustdirection --maxiterate 10 --add $prj/4_filter/"$virus"_draft_genomes."$run".fa $db/fasta/"$virus"_refs.fa > $prj/4_filter/"$virus"_draft_genomes."$run".ref_aligned.fa

#unalign draft genome alignment
sed '/^>/! s/\-//g' $prj/4_filter/"$virus"_draft_genomes."$run".ref_aligned.fa | sed 's/_R_//g' | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' > $prj/4_filter/"$virus"_draft_genomes."$run".ref_unaligned.fa

#retrieve draft denovo genomes
seqtk subseq $prj/4_filter/"$virus"_draft_genomes."$run".ref_unaligned.fa $prj/4_filter/"$virus"_draft_genomes."$run".list > $prj/4_filter/"$virus"_draft_genomes."$run".ref_unaligned.reoriented.fa

#create folder for blasting
mkdir $prj/5_remap

#retrieve individual draft sequences for mapping
while read i; do
grep -A1 ">Draft_${i}_" $prj/4_filter/"$virus"_draft_genomes."$run".ref_unaligned.reoriented.fa > $prj/5_remap/"$i".draft.fa
done < $prj/IDs.list

#map trimmed reads to draft genome
while read i; do
bbmap.sh -Xmx"$m"g threads="$t" maxindel=200 minid=0.98 in=$prj/1_reads/"$i"_trim_1.fastq.gz in2=$prj/1_reads/"$i"_trim_2.fastq.gz outm=$prj/5_remap/"$i".remapped.bam ref=$prj/5_remap/"$i".draft.fa
samtools sort -@ "$t" -o $prj/5_remap/"$i".remapped.sorted.bam $prj/5_remap/"$i".remapped.bam
rm $prj/5_remap/"$i".remapped.bam

done < $prj/IDs.list

