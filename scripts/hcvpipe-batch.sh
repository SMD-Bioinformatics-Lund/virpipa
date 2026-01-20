#!/bin/bash

subsample=1000000
outdir=/fs1/jonas/hcv/results/
logdir=/fs1/jonas/hcv/logs/
partition=low
dry=''

showhelp(){
    cat << EOF

Usage: $(basename "$0") [OPTIONS] -i <inputdir> | -c <csvfile>

A batch runner for the HCV pipeline that submits jobs via SLURM (sbatch).
It can either scan a directory for R1 FastQ files or use a CSV file for mapping.

REQUIRED ARGUMENTS (Choose one):
  -i, --inputdir PATH    Directory containing FastQ files (*R1*gz).
  -c, --csv PATH         CSV file with columns: clarity_sample_id, read1, sample_name.
                         (Takes precedence over --inputdir)

GENERAL OPTIONS:
  -o, --outdir PATH      Output directory (Default: $outdir)
  -l, --logdir PATH      Log directory (Default: $logdir)
  -p, --partition STR    SLURM partition (Default: $partition)
  -s, --subsample INT    Number of reads to subsample (Default: $subsample)

TOOLS:
  -n, --dryrun           Print the sbatch commands without executing them.
  -d, --debug            Enable verbose output and dump CSV contents.
  -h, --help             Display this help message and exit.

EXAMPLES:
  # Run using a directory of FastQ files:
  $(basename "$0") -i /path/to/fastqs -p high

  # Run using a CSV file (Dry run):
  $(basename "$0") -c samples.csv --dryrun

EOF
}

if [[ $# -eq 0 ]] ; then
	showhelp
	exit
fi

readopts=$(getopt -o hndo:s:l:p:i:c: --long help,dryrun,debug,outdir:,subsample:,logdir:,partition:,input:,csv: -n 'error' -- "$@")
eval set -- "$readopts"
dry=''

while true ; do
	case "$1" in
		-h|--help)
			showhelp
			exit 1 ;;
		-n|--dryrun)
			dry='echo'
			shift ;;
		-d|--debug)
			debug='1'
			shift ;;
		-o|--outdir)
			outdir="$2"
			shift 2 ;;
		-s|--subsample)
			subsample="$2"
			shift 2 ;;
		-l|--logdir)
			logdir="$2"
			shift 2 ;;
		-p|--partition)
			partition="$2"
			shift 2 ;;
		-i|--inputdir)
			inputdir="$2"
			shift 2 ;;
		-c|--csv)
			csv="$2"
			shift 2 ;;
		--) shift ; break ;;
		*) echo "unknown parameter" ; exit 1 ;;
	esac
done


sbatch_sample(){
	local jobname="$1"
	local sample="$2"
	local lid="$3"
	if [[ ! -z "$lid" ]] ; then
		$dry sbatch -J HCV-${jobname} --partition $partition $(dirname $0)/hcvpipe.sh -s $subsample -l $lid -o ${outdir}/${runname} $sample
	else
		$dry sbatch -J HCV-${jobname} --partition $partition $(dirname $0)/hcvpipe.sh -s $subsample -o ${outdir}/${runname} $sample
	fi	
}

runraw(){
	# Uses a directory with fastq files
	fulldir=$(readlink -f $inputdir)
	runname=$(basename $fulldir)
	for sample in ${inputdir}/*R1*gz ; do
		jobname=$(basename $sample)
		jobname=${jobname%%_*}
		if [[ -f "$fulldir"/../"$runname".tsv ]] ; then
			lid=$(grep $jobname "$fulldir"/../"$runname".tsv | cut -f2)
		fi
		sbatch_sample "$jobname" "$sample" "$lid"
	done
}

runcsv(){
	# Uses a csv file, parses the header and then created the sbatch commands
	declare -A csvdata
	exec 3< "$csv"
	IFS=',' read -r -a keys <&3
	while IFS=',' read -r -a values <&3 ; do
		#skip empty
		[[ -z "${values[0]}" ]] && continue
		for i in "${!keys[@]}"; do
			# Remove potential carriage returns (\r) from Windows-style CSVs
			key=$(echo "${keys[$i]}" | tr -d '\r')
			val=$(echo "${values[$i]}" | tr -d '\r')
			csvdata["$key"]="$val"
		done
		if [[ ! "$debug" == '' ]] ; then
			echo CSV contents:
			for k in "${!csvdata[@]}"; do
				echo "$k: ${csvdata[$k]}"
			done
			echo
		fi
		# send the paramteres
		sbatch_sample ${csvdata["clarity_sample_id"]} ${csvdata["read1"]} ${csvdata["sample_name"]}
	done
}

### main program
if [[ ! -d $outdir ]] ; then
	echo $outdir is not a directory
	exit
elif [[ -z "$inputdir" ]] && [[ -z "$csv" ]] ; then
	echo You must provide either inputdir or csv
	exit
elif [[ ! -z "$inputdir" ]] && [[ ! -z "$csv" ]] ; then
	echo You cannot provide both inputdir and csv
	exit
fi

cd $logdir

if [[ ! -z "$csv" ]] ; then
	runcsv
else
	runraw
fi

cd - > /dev/null
