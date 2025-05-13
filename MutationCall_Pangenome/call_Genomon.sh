#! /bin/sh

blatDIR=/home/package/blat/linux.x86_64/
export PATH=${blatDIR};${PATH}
samtoolsDIR=/usr/local/package/samtools/1.9/bin/
export PATH=${samtoolsDIR};${PATH}
fasta=/home/reference/Human_Pangenome_Reference_Consortium/Graphs/GRCh38/AF_filtered/hprc-v1.1-mc-grch38.d9.gbz.GRCh38referencepaths.fa
virtualenvGenomon=/home/usr/Genomon2/bin
source ${virtualenvGenomon}

#Aux Function
is_file_exists() {
	if [ -f $1 ]; then
		echo "$1 exists."
		return 0
	fi
	echo "$1 does not exists."
	exit 1
}

check_mkdir() {
	if [ -d $1 ]; then
		echo "$1 exists."
	else
		echo "$1 does not exits."
		mkdir -p $1
	fi
}


#Genomon Config
base_quality=15
min_allele_freq=0.005
max_allele_freq=0.10
min_depth=8
min_variant_read=5
fisher_value=1.0
samtools_params="-q 0 -BQ0 -d 10000000 --ff UNMAP,QCFAIL"
cutNnumFisher=4
score_difference=5
window_size=200
max_depth=5000
exclude_sam_flags=3328
EBcall_mapq=0
EBcall_baseq=15


usage() {
    echo "Arguments=[information table] [-f | --fasta [ARG]] [(-d | --directory  [ARG])] [(-p | --project [ARG])] [(-T | --tag [ARG])][(-o | --output [ARG])][(-j | --jobfile [ARG])][(-J)]"
    echo "Options:"
    echo "  -h, --help"
    echo "  -f | --fasta [ARG] (mandatory), reference fasta file, default = ${fasta}"
    echo "  -d | --directory [ARG] (option), output parent directory, default = ${DIR}"
    echo "  -p | --project [ARG] (option), project names, default = date"
    echo "  -T | --tag [ARG] (option), output TAG, default = file name"
    echo "  -o | --output [ARG] (option), direct path of output. if -o is used, -T, -p, and -d are ignored, default = NULL"
    echo "  -j | --jobfile [ARG] (option), job list, default = ${jobfile}"
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
		'-f'|'--fasta' )
			if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
				echo "option requires an argument $1" 1>&2
				exit 1
			fi
			fasta=$2
			shift 2
			;;
		'-c' )
			if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
				echo "option requires an argument $1" 1>&2
				exit 1
			fi
			genomonconfig=$2
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

TBAM=${param[0]} ; is_file_exists ${TBAM}
NBAM=${param[1]} ; is_file_exists ${NBAM}
OUTPUTTAG=${param[2]} ; check_mkdir $(dirname ${OUTPUTTAG})
Panel=${param[3]} ; is_file_exists ${Panel}

##check files
is_file_exists ${genomonconfig}

##Fisher call
fisher comparison -1 ${TBAM} -2 ${NBAM} -o ${OUTPUTTAG}.U1snRNA.Genomon.fisher.txt --ref_fa ${fasta} --samtools_path $(which samtools) --base_quality ${base_quality} --min_allele_freq ${min_allele_freq} --max_allele_freq ${max_allele_freq} --min_depth ${min_depth} --min_variant_read ${min_variant_read} --fisher_value ${fisher_value} --samtools_params "${samtools_params}" || exit $?

##Local realignment by blat
mutfilter realignment -s ${score_difference} -w ${window_size} -t ${OUTPUTTAG}.U1snRNA.Genomon.fisher.txt -1 ${TBAM} -2 ${NBAM} -o ${OUTPUTTAG}.U1snRNA.Genomon.realign.txt -r ${fasta} -b $(which blat) --exclude_sam_flags ${exclude_sam_flags} || exit $?
rm ${OUTPUTTAG}.U1snRNA.Genomon.fisher.txt|| exit $?

#Empirical Baysian mutation calling

EBFilter -f anno -q ${EBcall_mapq} -Q ${EBcall_baseq} --ff 4 ${OUTPUTTAG}.U1snRNA.Genomon.realign.txt ${TBAM} ${Panel} ${OUTPUTTAG}.U1snRNA.Genomon.EB.tsv --loption || exit $?
rm ${OUTPUTTAG}.U1snRNA.Genomon.realign.txt || exit $?
deactivate
