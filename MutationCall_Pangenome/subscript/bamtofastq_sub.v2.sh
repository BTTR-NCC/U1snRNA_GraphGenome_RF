#! /bin/bash

#Aux Function
is_file_exists() {
	if [ -f $1 ]; then
		echo "$1 exists."
		return 0
	fi
	echo "$1 does not exists."
	exit 1
}

bbDIR=/home/package/biobambam/2.0.146/bin/
samtoolsDIR=/usr/local/package/samtools/1.9/bin/
export PATH=${bbDIR};${PATH}
export PATH=${samtoolsDIR};${PATH}

rminput="FALSE"
rmunpair="FALSE"
cram="FALSE"

usage() {
    echo "Arguments=[BAM] [(-o | --output [ARG])] [(-j | --jobfile [ARG])] [(-J)]"
    echo "Options:"
    echo "  -h, --help"
    echo "  -c | --cram [ARG] (option), set when the input file is cram. default = FALSE"
    echo "  -f | --fasta [ARG] (option), reference fasta file, must be needed when the input file is cram. default = FALSE"
    echo "  -o | --output [ARG] (option), output path. default = BAM file name"
    echo "  -j | --jobfile [ARG] (option), job list, default = ${JOBLIST}"
    echo "  -J (option), do not make a job list, default = FALSE"
	echo
    exit 1
}
for OPT in "$@"
do
	case "$OPT" in
		'-h'|'--help' )
			usage
			exit 1
			;;
		'--rminput' | '--replace' )
			rminput="TRUE"
			shift 1
			;;
		'--paironly' | '--rmunpair' )
			rmunpair="TRUE"
			shift 1
			;;
		'-c' | '--cram' )
			cram="TRUE"
			shift 1
			;;
		'-f'|'--fasta' )
			if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
				echo "option requires an argument $1" 1>&2
				exit 1
			fi
			fasta=$2
			shift 2
			;;
		'-' )
			shift 1
			param+=( "$@" )
			break
			;;
		-*)
			echo "illegal option -- '$(echo $1 | sed 's/^-*//')'" 1>&2
			exit 1
			;;
		*)
			if [[ ! -z "$1" ]] && [[ ! "$1" =~ ^-+ ]]; then
				#param=( ${param[@]} "$1" )
				param+=( "$1" )
				shift 1
			fi
			;;
	esac
done

#input data
BAM=${param[0]}; is_file_exists ${BAM}
OUTPUT=${param[1]};

samtools quickcheck ${BAM} || exit 1

inputT=${OUTPUT}_temp.txt
inputS=${OUTPUT}_single_end_output.txt.gz
inputO=${OUTPUT}_unmatched_first_output.txt.gz
inputO2=${OUTPUT}_unmatched_second_output.txt.gz
inputF=${OUTPUT}_sequence1.fq.gz
inputF2=${OUTPUT}_sequence2.fq.gz

if [[ "${cram}" == "TRUE" ]]; then
  if [[ -f "${fasta}" ]]; then
    bamtofastq gz=1 inputformat=cram reference=${fasta} collate=1 exclude=SECONDARY,SUPPLEMENTARY T=${inputT} S=${inputS} O=${inputO} O2=${inputO2} filename=${BAM} F=${inputF} F2=${inputF2} || exit $?
  else
    echo "Reference fasta file must be specified when using cram input."
    exit 1
  fi
else
  bamtofastq gz=1 collate=1 exclude=SECONDARY,SUPPLEMENTARY T=${inputT} S=${inputS} O=${inputO} O2=${inputO2} filename=${BAM} F=${inputF} F2=${inputF2} || exit $?
fi

if [[ -f "${inputF}" ]] && [[ -f "${inputF2}" ]] ; then
	if [[ "${rmunpair}" == "TRUE" ]] ; then
		rm ${inputS} ${inputO} ${inputO2}
	fi
	if [[ "${rminput}" == "TRUE" ]] ; then
		rm ${BAM}
	fi
fi
