#!/usr/bin/env bash
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=16
#SBATCH --partition=grace-normal

ml apptainer/1.0.2 GCC/8.3.0 seqtk/1.3 pigz/2.4 datamash/1.5 Python/3.7.4

# define defaults
cpus=16
type=r11b2L25
scripts="/fs1/jonas/hcv/scripts"
containerdir='/fs1/resources/containers'
refdir='/fs1/jonas/hcv/refgenomes'
genome_files=( "${refdir}/1a-AF009606.fa" "${refdir}/1a-M62321.fa" "${refdir}/1b-D90208.fa" "${refdir}/2a-D00944.fa" "${refdir}/2b-D10988.fa" "${refdir}/3a-D17763.fa" "${refdir}/3k-HPCJK049E1.fa" "${refdir}/4a-GU814265.fa" "${refdir}/5a-Y13184.fa" "${refdir}/6a-Y12083.fa" "${refdir}/6g-HPCJK046E2.fa" "${refdir}/7a-EF108306.fa" )
# genome_files=( "${refdir}/1a-AF009606.fa" "${refdir}/1b-D90208.fa" )

# sent="apptainer exec -B /fs1,/local ${containerdir}/sentieon_202308.01.sif sentieon"
#sent="apptainer exec -B /fs1,/local ${containerdir}/sentieon_202503--7e7ce56c0d0199c5.sif"
sent="apptainer exec -B /fs1,/local ${containerdir}/sentieon_202308.03.sif sentieon"
samt="apptainer exec -B /fs1,/local ${containerdir}/samtools_1.21.sif samtools"
bgz="apptainer exec -B /fs1,/local ${containerdir}/samtools_1.21.sif bgzip"
bcft="apptainer exec -B /fs1,/local ${containerdir}/bcftools_1.21.sif bcftools"
mum="apptainer exec -B /fs1,/local ${containerdir}/mummer3.23.sif"
mafftbin="apptainer exec -B /fs1,/local ${containerdir}/mafft7.525.sif mafft"
spadesbin="apptainer exec -B /fs1,/local ${containerdir}/spades_3.15.5.sif spades.py"
hostilebin="apptainer exec -B /fs1,/local ${containerdir}/hostile_1.1.0.sif hostile"
freebayes="apptainer exec -B /fs1,/local ${containerdir}/freebayes_1.3.8.sif freebayes"
blastn="apptainer exec -B /fs1,/local ${containerdir}/blast_2.16.0.sif blastn"
pyplot="apptainer exec -B /fs1,/local ${containerdir}/pyplot.sif python ${scripts}/kderug.py"

export LC_NUMERIC=en_US.UTF-8    # otherwise datamash recognizes , as decimal separator
export HOSTILE_CACHE_DIR=/fs1/resources/ref/micro/hostile

#echo $id
#echo $genome_files

function showshorthelp() {
	echo "hcvpipe.sh -h/--help for help"
	echo
}

function showhelp() {
	echo
	echo 'Syntax: hcvpipe.sh [OPTIONS] read1 [read2]'
	echo
	echo 'Simple HCV pipeline with optional subsampling'
	echo
	echo 'Positional arguments:'
	echo '    read1 <read2>             Give at least a forward read as an argument'
	echo '                              if only forward is given, reverse will be inferred'
	echo '                              but you probably want to change this'
	echo 'Optional arguments:'
	echo '   -s <i>, --subsample <n>    Subsample reads. [default 500 000]'
	echo '   -o <s>, --outdir <s>       Sets the root dir for output folders'
	echo '   -r <s>, --reference <s>    Use this reference'
	echo '   -c <i>, --cpus <i>         Number of CPUs to use. [default 16]'
	echo '   -H, --hostile              Do NOT remove human reads with hostile'
	echo '   -n, --dry-run              Dry run'
	echo '   -f, --force                Force rerun'
	echo '   -h, --help                 This help'
	echo
}

# read the options

if [[ -z "$1" ]]; then
	echo "You must provide arguments."
	showshorthelp
	exit
fi


readopts=$(getopt -o hnfc:s:Ho:r: --long help,dryrun,force,cpus,subsample:,hostile,outdir:,reference: -n 'error' -- "$@")
#echo $readopts
eval set -- "$readopts"
pc=''
force=0
outdirroot=""
hostile='true'

while true ; do
	case "$1" in
		-s|--subsample)
			case "$2" in
				"") subsamplereads='500000' ;shift 2 ;;
				*) subsamplereads=$2 ;shift 2 ;;
			esac ;;
		-o|--outdir)
			outdirroot="$2"
			shift 2;;
		-c|--cpus)
			cpus="$2"
			shift 2;;
		-r|--reference)
			genome_files=( "$2" )
			shift 2;;
		-n|--dryrun)
			pc=echo
			shift;;
		-H|--dryrun)
			hostile='false'
			shift;;
		-f|--force)
			force=1
			shift;;
		-h|--help) showhelp; exit 0;;
		--) shift ; break ;;
		*) echo "Unknown parameter!" ; exit 1 ;;
	esac
done

echo precommand is $pc
echo subsample to $subsamplereads

# parse the positional arguments

if [[ $# -gt 2 ]] ; then
	echo 'Wrong number of arguments'
	exit
elif [[ $# -eq -0 ]] ; then
	echo 'You must supply at least a forward read'
	exit
elif [[ $# -eq 1 ]] ; then
	# parse R1 and infer R2
	r1org="$1"
	r2org="${r1org//_R1_/_R2_}"
else
	r1org="$1"
	r2org="$2"
fi

r1base=$(basename $r1org)
r2base=$(basename $r2org)
id=${r1base%%_*}

# set name of output dir in rootdir

if [[ $outdirroot == "" ]] ; then
	outdir=${id}
else
	outdir="${outdirroot%/}/${id}"
fi

### FUNCTIONS ###

function abspath() {
	if [[ -d $1 ]]; then
		cd $1
		printf "%s\n" "$(pwd)"
	else
		cd "$(dirname "$1")"
		printf "%s/%s\n" "$(pwd)" "$(basename "$1")"
	fi
}

function subsample() {
	local read1="$1"
	local read2="$2"
	local nofreads=$3
	local outsubdir=$4
	echo precommand is $pc
	if [[ $pc == 'echo' ]] ; then
		echo subsampling
		echo seqtk sample -s100 $read1 $nofreads \| pigz \> ${outdir}/${outsubdir}/${r1base/.fastq/.sub.fastq}
		echo seqtk sample -s100 $read2 $nofreads \| pigz \> ${outdir}/${outsubdir}/${r2base/.fastq/.sub.fastq}
	else
		if [[ ! -f ${outdir}/fastq/${r1base/.fastq/.sub.fastq} ]] ; then
			seqtk sample -s100 $read1 $nofreads | pigz > ${outdir}/${outsubdir}/${r1base/.fastq/.sub.fastq}
		fi
		if [[ ! -f ${outdir}/fastq/${r2base/.fastq/.sub.fastq} ]] ; then
			seqtk sample -s100 $read2 $nofreads | pigz > ${outdir}/${outsubdir}/${r2base/.fastq/.sub.fastq}
		fi
	fi
}

function umimap() { # create bam files while first making a consensus using the UMIs
	local localref="$1"
	local shortrefname="$2"

	bamfile="${outdir}/bam/${id}-${shortrefname}.${type}.bwa.umi.sort.bam"
	echo r1 $r1
	echo r2 $r2

	# consensus from UMI

	$sent umi extract -d 3M2S+T,3M2S+T $r1 $r2 \
		| $sent bwa mem \
		-R "@RG\\tID:$id\\tSM:$id\\tLB:$id\\tPL:illumina" \
		-t ${cpus} \
		-k 11 \
		-B 2 \
		-L 25 \
		-p -C $localref - \
		| $sent umi consensus --copy_tags XR,RX,MI,XZ -o ${outdir}/tmp/consensus.fastq.gz

	#		pipe to tee -a noumi.sam    # if you wish to create the non-consensus bam

	$sent bwa mem \
		-R "@RG\\tID:$id\\tSM:$id\\tLB:$id\\tPL:illumina" \
		-t ${cpus} \
		-k 11 \
		-B 2 \
		-L 25 \
		-p -C $localref ${outdir}/tmp/consensus.fastq.gz \
		| $sent util sort -i - \
		-o ${bamfile} \
		--sam2bam --umi_post_process

	#	$sent util sort -i noumi.sam -o ${id}-${genome_file_short}.${type}.bwa.sort.bam --sam2bam	# for non-consensus bam
	#	rm noumi.sam

	filterbam=${bamfile/.sort.bam/.filter.sort.bam}
	$samt view -@ 30 ${bamfile} -e "sclen < 30" --with-header --bam --output ${filterbam}
	$samt index ${bamfile}
	$samt index ${filterbam}
	$samt stats ${filterbam} > ${filterbam}.stats
}

function bam2fasta() { # gets a fasta file from bam with bcftools
	local localbam="$1"
	local localref="$2"
	local shortrefname="$3"
	local idref="${id}-${shortrefname}"
	$pc $bcft mpileup -Ob -d 10000 -o ${outdir}/bcf/${idref}.bcf -f ${localref} ${localbam}
	$bcft call --ploidy 2 --variants-only -Am -O z -o ${outdir}/vcf/${idref}.vcf.gz ${outdir}/bcf/${idref}.bcf
	$bcft index  ${outdir}/vcf/${idref}.vcf.gz
	$bcft consensus -f ${localref} ${outdir}/vcf/${idref}.vcf.gz >  ${outdir}/fasta/${idref}.fasta
	sed -i "s/>.*/>${idref}/" ${outdir}/fasta/${idref}.fasta
	$bcft stats ${outdir}/vcf/${idref}.vcf.gz > ${outdir}/vcf/${idref}.vcf.gz.stats
}

function getbestsubtype() {
	for statsfile in ${outdir}/bam/*stats ; do
		printf "${statsfile%%.*}\t" | sed 's:.*/::' | sed 's:-:\t:'
		grep error ${statsfile} | awk '{print $4}'
	done | \
		datamash -f -g1 min 3 | cut -f2 # get the lowest error from the generated bam.stats files
}

function createreport() {
	stats="$1"
	local id=${stats%.vcf.gz.stats}
	local id=$(echo $id | cut -d'-' -f1)
	ref=$(echo $id | cut -d'-' -f2,3)
	snps=$(sed -n 's/^SN\t0\tnumber of SNPs:\t\([0-9]*\)$/\1/p' $stats)
	multiallelic=$(sed -n 's/^SN\t0\tnumber of multiallelic SNP sites:\t\([0-9]*\)$/\1/p' $stats)
	af0=$(sed -n 's/^AF\t0\t0.000000\t\([0-9]*\)\t.*/\1/p' $stats)
	af99=$(sed -n 's/^AF\t0\t0.990000\t\([0-9]*\)\t.*/\1/p' $stats)
	echo '# VCF stats'
	printf "subtype\t${subtype%-*}\n"
	printf "reference\t${subtype}\n"
	printf "id\t${id}\n"
	printf "snps\t${snps}\n"
	printf "multiallelic_snps\t${multiallelic}\n"
	printf "AF0\t${af0}\n"
	printf "AF0.99\t${af99}\n"
}

### Main program

# Create outdirs

if [[ -d ${outdir} ]] && [[ $force -eq 0 ]] ; then
	echo Directory $outdir does already exist. force with -f
	exit
else
	$pc mkdir -p $outdir 2> /dev/null
	for i in fastq bam bcf vcf fasta spades mummer tmp results ; do
		$pc mkdir ${outdir}/${i} 2> /dev/null
	done
fi

# remove human reads
if [[ $hostile == 'true' ]] ; then
	$pc $hostilebin clean --offline --fastq1 $r1org --fastq2 $r2org --out-dir ${outdir}/fastq > ${outdir}/results/hostile.json
	r1=${outdir}/fastq/${r1base/.fastq.gz/.clean_1.fastq.gz}
	r2=${outdir}/fastq/${r2base/.fastq.gz/.clean_2.fastq.gz}
else
	r1=$r1org
	r2=$r2org
fi

# ugly hack since I stupidly used global r1 and r2
r1nosub=$r1
r2nosub=$r2


# non-optional subsampling to 250k reads for finding correct subtype
subsample $r1 $r2 250000 tmp
r1=${outdir}/tmp/${r1base/.fastq/.sub.fastq}
r2=${outdir}/tmp/${r2base/.fastq/.sub.fastq}


#loop through reference genomes and create bams and fasta
for genome_file in "${genome_files[@]}" ; do
	genome_file_short=$(basename $genome_file)
	genome_file_short=${genome_file_short%.fa}
	echo $id $genome_file_short
	if [[ $pc == 'echo' ]] || [[ -f ${outdir}/bam/${id}-${genome_file_short}.${type}.bwa.umi.sort.bam ]]; then
		echo umimap  "${genome_file}" "${genome_file_short}"
		echo bamtofasta "${outdir}/bam/${id}-${genome_file_short}.${type}.bwa.umi.filter.sort.bam" "${genome_file}" "${genome_file_short}"
		continue
	else
		[[ ! -f ${genome_file}.bwt ]] && $sent bwa index $genome_file
		umimap "${genome_file}" "${genome_file_short}"
		bam2fasta  "${outdir}/bam/${id}-${genome_file_short}.${type}.bwa.umi.filter.sort.bam" "${genome_file}" "${genome_file_short}"
	fi
done

# ugly restoration of non-subsampled reads
r1=$r1nosub
r2=$r2nosub

# optional subsampling
if [[ ! $subsamplereads == '' ]] ; then
	subsample $r1 $r2 ${subsamplereads} fastq
	r1=${outdir}/fastq/${r1base/.fastq/.sub.fastq}
	r2=${outdir}/fastq/${r2base/.fastq/.sub.fastq}
fi

# get the best subtype
subtype=$(getbestsubtype)

# copying results data and create cram from the best results
# mapping to the best subtype with all reads
$pc cp ${outdir}/fasta/${id}-${subtype}.fasta ${outdir}/results/
$pc $samt view -O cram,embed_ref -T ${refdir}/${subtype}.fa ${outdir}/bam/${id}-${subtype}.${type}.bwa.umi.filter.sort.bam -o ${outdir}/results/${id}-${subtype}.cram
$pc $samt index ${outdir}/results/${id}-${subtype}.cram
$pc cp ${outdir}/vcf/${id}-${subtype}.vcf.gz* ${outdir}/results

# de novo assembly (possibly on subsampled reads)
if [[ ! -f ${outdir}/spades/${id}.spades ]] ; then
	$pc $spadesbin -1 $r1 -2 $r2 --rnaviral -o ${outdir}/spades --threads ${cpus}
	$pc cp ${outdir}/spades/contigs.fasta ${outdir}/spades/${id}.spades
fi


# mapping back to the best hit - mainly for compensating for gaps in the mapping to the ref
$pc cd ${outdir}/tmp
$pc cp ${outdir}/results/${id}-${subtype}.fasta 4mafft.fa
$pc sed -i "s/>${id}-${subtype}/>${id}/" 4mafft.fa

mkdir fastasplit
$pc cd fastasplit
$pc awk '/^>/ {OUT=substr($0,2) ".fa"}; OUT {print >OUT}' ${outdir}/spades/${id}.spades
$pc cd -
$pc $mum nucmer -p ${id} ${outdir}/fasta/${id}-${subtype}.fasta ${outdir}/spades/${id}.spades
$pc $mum show-tiling ${id}.delta | sed -n '2,$p' | cut -f7,8 > ${outdir}/tmp/${id}.tiling
while read -r direction contig ; do
	if [[ $direction == '-' ]] ; then
		seqtk seq -r fastasplit/${contig}.fa >> 4mafft.fa
	else
		cat fastasplit/${contig}.fa >> 4mafft.fa
	fi
done < ${outdir}/tmp/${id}.tiling
# align contigs and do not use iupac
$pc $mafftbin 4mafft.fa > ${id}.mafft 
python $scripts/consensus_fasta.py ${id}.mafft 2> ${outdir}/mummer/${id}-denovocons.out > ${outdir}/mummer/${id}-denovocons.fasta
python $scripts/consensus_fasta.py ${id}.mafft gapped 2> /dev/null > ${outdir}/mummer/${id}-denovocons-gapped.fasta
$pc $samt faidx ${outdir}/mummer/${id}-denovocons-gapped.fasta
$pc $sent bwa index "${outdir}/mummer/${id}-denovocons-gapped.fasta"
$pc $samt faidx ${outdir}/mummer/${id}-denovocons.fasta
$pc $sent bwa index "${outdir}/mummer/${id}-denovocons.fasta"
$pc cd -

# convert IUPC to random nuc in order to work with freebayes
seqtk randbase ${outdir}/mummer/${id}-denovocons.fasta > ${outdir}/mummer/${id}-denovocons-randbase.fasta
$sent bwa index ${outdir}/mummer/${id}-denovocons-randbase.fasta
# map to new fasta starting from the gapped fasta from aligned de novo contigs. This is to generate a new fasta template
$pc umimap "${outdir}/mummer/${id}-denovocons-randbase.fasta" "denovocons-randbase"
# $pc umimap "${outdir}/mummer/${id}-denovocons.fasta" "denovocons-iter1"
# maybe do this with ploidy 1 instead of 2 to not have ambiguous calls? Plus one with IUPAC codes for say 10%?
# 250417 change the min-alternative-cutoff to 15% since this seems to be the global consensus
$freebayes --ploidy 1 --min-coverage 3 --min-base-quality 20 --min-alternate-fraction 0.15 --min-mapping-quality 60 -f ${outdir}/mummer/${id}-denovocons-randbase.fasta ${outdir}/bam/${id}-denovocons-randbase.${type}.bwa.umi.filter.sort.bam > ${outdir}/vcf/${id}-freebayes-randbase.vcf
$freebayes --ploidy 2 --min-coverage 3 --min-base-quality 20 --min-alternate-fraction 0.15 --min-mapping-quality 60 -f ${outdir}/mummer/${id}-denovocons-randbase.fasta ${outdir}/bam/${id}-denovocons-randbase.${type}.bwa.umi.filter.sort.bam > ${outdir}/vcf/${id}-freebayes-randbase-iupac.vcf
# index the vcf
$bgz -f -i ${outdir}/vcf/${id}-freebayes-randbase.vcf
$bcft index ${outdir}/vcf/${id}-freebayes-randbase.vcf.gz
$bgz -f -i ${outdir}/vcf/${id}-freebayes-randbase-iupac.vcf
$bcft index ${outdir}/vcf/${id}-freebayes-randbase-iupac.vcf.gz
# get the new fasta based on freebayes calling
#$samt index ${outdir}/mummer/${id}-denovocons-randbase.fasta
$bcft consensus -f ${outdir}/mummer/${id}-denovocons-randbase.fasta ${outdir}/vcf/${id}-freebayes-randbase.vcf.gz > ${outdir}/fasta/${id}-freebayes.fasta
$bcft consensus -f ${outdir}/mummer/${id}-denovocons-randbase.fasta ${outdir}/vcf/${id}-freebayes-randbase-iupac.vcf.gz > ${outdir}/fasta/${id}-freebayes-iupac.fasta

# dirty copy to results and create report FIXME
cp ${outdir}/fasta/${id}-freebayes-iupac.fasta ${outdir}/results/
awk -vFS="" 'NR>1 {for(i=1;i<=NF;i++)w[toupper($i)]++}END{for(i in w) print i,w[i]}' ${outdir}/results/${id}-freebayes-iupac.fasta | sort -nr -k2 > ${outdir}/results/${id}-freebayes-iupac.fastanucfreq.tsv

# and then map once more to the file without IUPAC and vcf call and that will give you files to use
$samt faidx ${outdir}/fasta/${id}-freebayes.fasta
$sent bwa index ${outdir}/fasta/${id}-freebayes.fasta
$pc umimap "${outdir}/fasta/${id}-freebayes.fasta" "freebayes"

# create various tracks for igv
for minfrac in 0.01 0.05 0.1 0.15 0.2 0.3 0.4; do
	$freebayes --ploidy 2 --min-coverage 3 --min-base-quality 20 --min-alternate-fraction ${minfrac} --min-mapping-quality 60 -f ${outdir}/fasta/${id}-freebayes.fasta ${outdir}/bam/${id}-freebayes.${type}.bwa.umi.filter.sort.bam > ${outdir}/vcf/${id}-freebayes-m${minfrac}.vcf
	$bgz -f -i ${outdir}/vcf/${id}-freebayes-m${minfrac}.vcf
	$bcft index ${outdir}/vcf/${id}-freebayes-m${minfrac}.vcf.gz
	$bcft stats ${outdir}/vcf/${id}-freebayes-m${minfrac}.vcf.gz > ${outdir}/vcf/${id}-freebayes-m${minfrac}.vcf.gz.stats
done

# render cram for consensus mapped bam
echo three
$pc $samt view -O cram,embed_ref -T ${outdir}/fasta/${id}-freebayes.fasta ${outdir}/bam/${id}-freebayes.${type}.bwa.umi.filter.sort.bam -o ${outdir}/results/${id}-freebayes.cram
$pc $samt index ${outdir}/results/${id}-freebayes.cram
$pc cp ${outdir}/vcf/${id}-freebayes-m*.vcf* ${outdir}/results/
$pc cp ${outdir}/fasta/${id}-freebayes.fasta ${outdir}/results/
$pc cp ${outdir}/fasta/${id}-freebayes.fasta.fai ${outdir}/results/


# create report files with stats and quality FIXME! freebayes part
ln -s ${outdir}/vcf/${id}-freebayes-m0.15.vcf.gz.stats ${outdir}/vcf/${id}-freebayes.vcf.gz.stats
for report in ${id}-${subtype} ${id}-freebayes ; do
	createreport ${outdir}/vcf/${report}.vcf.gz.stats > ${outdir}/results/${report}.report.tsv
	echo '# COVERAGE' >> ${outdir}/results/${report}.report.tsv
	$samt coverage ${outdir}/results/${report}.cram | datamash transpose | sed -n '2,$p' >> ${outdir}/results/${report}.report.tsv
	awk -vFS="" 'NR>1 {for(i=1;i<=NF;i++)w[toupper($i)]++}END{for(i in w) print i,w[i]}' ${outdir}/results/${report}.fasta | sort -nr -k2 > ${outdir}/results/${report}.fastanucfreq.tsv
done

# blast the final fasta and the de novo
# cp ${outdir}/spades/${id}.spades ${outdir}/results/
for fasta in ${outdir}/results/${id}-freebayes.fasta  ${outdir}/spades/${id}.spades ; do
    printf "query acc.ver\tsubject acc.ver\t%% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. star\
t\ts. end\tevalue\tbit score\n" > ${fasta}.blast
    $blastn -query $fasta -db ${refdir}/hcvglue/hcvgluerefs -outfmt 6 >> ${fasta}.blast
done

# plot kde + rug of vcf data with 1% cutoff
cd ${outdir}/results/
$pc zcat ${outdir}/vcf/${id}-freebayes-m0.1.vcf | \
	grep 'GT:DP:AD:RO:QR:AO:QA:GL' | \
	sed 's/.*\t//' | \
	cut -d':' -f2,4 | \
	tr ':' '\t' | \
	awk '{print $2/$1}' > ${outdir}/tmp/${id}.mixin
$pc $pyplot ${outdir}/tmp/${id}.mixin
cd -

# log genome coverage at 1x 10x 100x 1000x
printf "id\t1x\t10x\t100x\t1000x\n${id}\t" > ${outdir}/results/${id}-coverage.tsv

$samt depth -r ${id}:100-9600 ${outdir}/results/${id}-freebayes.cram | \
awk 'BEGIN { total=9501; cov1=0; cov10=0; cov100=0; cov1000=0 }
     { if ($3 >= 1) cov1++; if ($3 >= 10) cov10++; if ($3 >= 100) cov100++; if ($3 >= 1000) cov1000++ }
     END {
         printf "%.2f\t", (cov1/total)*100;
         printf "%.2f\t", (cov10/total)*100;
         printf "%.2f\t", (cov100/total)*100;
         printf "%.2f\n", (cov1000/total)*100;
     }' >> ${outdir}/results/${id}-coverage.tsv
