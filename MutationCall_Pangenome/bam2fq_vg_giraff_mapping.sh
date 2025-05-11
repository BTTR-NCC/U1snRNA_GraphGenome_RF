#! /bin/sh
#$ -S /bin/bash
#$ -e /home/ha6434/command/vg/log/
#$ -o /home/ha6434/command/vg/log/
#$ -cwd
#$ -l s_vmem=24G,mem_req=24G
#$ -pe def_slot 6

source /home/ha6434/util/utility.sh
source /home/ha6434/util/source_share.sh

module load vg/1.53.0 samtools/1.9

bamtofastq=/home/ha6434/command/biobambam/subscript/bamtofastq_sub.v2.sh
rmfastq="FALSE"
mkbam="TRUE"
maxlen=3000
regionbed=""
inputbam=""

#clone of bam2fq_vg_giraff_mapping_sub.v2.3.sh for github

#parameter
compress="FALSE"

usage() {
	echo "Arguments=[tumor BAM] [tumor TAG] [normal bam] [(-d | --divconfig  [ARG])] [(-o | --output [ARG])] [(-c | --genomonconfig  [ARG])] [(-j | --jobfile [ARG])][(-J)]"
	echo "Options:"
	echo "  -h, --help"
	echo "  -B | --nobam [ARG] (option), do not produce bam, default = ${mkbam}"
	echo "  -b | --bed [ARG] (option), bam2fastq for selected bam region, default = """
	echo "  -r | --reference [ARG] (option), reference file for mapping with a suffix of .gbz, default = None"
	echo "  -m | --maxlen [ARG] (option), reads with fragment lengths greater than N will not be marked properly paired in SAM/BAM/CRAM, default = ${maxlen}"
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
		'-B'|'--nobam' )
			mkbam="FALSE"
			shift 1
			;;
		'-b'|'--bed' )
			if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
				echo "option requires an argument $1" 1>&2
				exit 1
			fi
			regionbed=$2
			shift 2
			;;
		'-r'|'--reference' )
			if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
				echo "option requires an argument $1" 1>&2
				exit 1
			fi
			reference=$2
			shift 2
			;;
		'-m'|'--maxlen' )
			if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
				echo "option requires an argument $1" 1>&2
				exit 1
			fi
			maxlen=$2
			shift 2
			;;
		'-f'|'--fasta' )
			if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
				echo "option requires an argument $1" 1>&2
				exit 1
			fi
			fasta=$2
			shift 2
			;;
		'-x'|'--xg' )
			if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
				echo "option requires an argument $1" 1>&2
				exit 1
			fi
			xg=$2
			shift 2
			;;
		'--rm' )
			rmfastq="TRUE"
			shift 1
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

BAM=${param[0]} ; is_file_exists ${BAM}
OUTPUTTAG=${param[1]}

if [[ -z "${OUTPUTTAG}" ]] ; then
	echo "OUTPUT (Argument 2) is not specified."
	exit 1
else
	check_mkdir $(dirname ${OUTPUTTAG})
fi

if [[ -z "${reference}" ]]; then
	"Reference graph for mapping (.gbz or .gg) is necessary. Please use -r for specity."
	exit 1
else
	echo "Check reference files..."
	is_file_exists ${reference}
	minimizer=$(echo ${reference} | sed -e  "s/.giraffe.gbz\|.gbz\|.gg//g").min ; is_file_exists ${minimizer}
	dist=$(echo ${reference} | sed -e  "s/.giraffe.gbz\|.gbz\|.gg//g").dist ; is_file_exists ${dist}
	if [[ ${reference} =~ .gg$ ]]; then
		gbwt=$(echo ${reference} | sed -e  "s/.gg//g").gbwt ; is_file_exists ${gbwt}
		xg=$(echo ${reference} | sed -e  "s/.gg//g").xg ; is_file_exists ${xg}
	fi
fi

#bamtofastq
if [[ ${BAM} =~ bam$ ]] ; then
	if [[ -n "${regionbed}" ]]; then
		inputbam=$(echo ${BAM} | sed -e "s/.bam\|.cram//g")_region.bam
		samtools view ${BAM} -L ${regionbed} -b -o ${inputbam} || exit $?
		sh ${bamtofastq} ${inputbam} ${OUTPUTTAG} || exit $?
	else
		sh ${bamtofastq} ${BAM} ${OUTPUTTAG} || exit $?
	fi
elif [[ ${BAM} =~ cram$ ]] ; then
	if [[ -z "${fasta}" ]]; then
		echo "Reference fasta for cram file must be specified when using cram input."
		exit 1
	else
		if [[ -n "${regionbed}" ]]; then
			inputbam=$(echo ${BAM} | sed -e "s/.bam\|.cram//g")_region.bam
			samtools view ${BAM} -T ${fasta} -L ${regionbed} -b -o ${inputbam} || exit $?
			sh ${bamtofastq} ${inputbam} ${OUTPUTTAG} || exit $?
		else
			sh ${bamtofastq} ${BAM} ${OUTPUTTAG} -c -f ${fasta}|| exit $?
		fi
	fi
else
  	echo "Invalid imput file. Check the imput file type is bam or cram."
  	exit 1
fi

fastq1=${OUTPUTTAG}_sequence1.fq.gz ; is_file_exists ${fastq1}
fastq2=${OUTPUTTAG}_sequence2.fq.gz ; is_file_exists ${fastq2}
rm ${OUTPUTTAG}_single_end_output.txt.gz ${OUTPUTTAG}_unmatched_first_output.txt.gz ${OUTPUTTAG}_unmatched_second_output.txt.gz
if [[ -n "${inputbam}" ]] && [[ -f "${inputbam}" ]]; then
	rm ${inputbam}
fi

#get RGinfo
TAG=$(basename $(dirname ${BAM}))
RG=$(samtools view -H ${BAM} | grep ^@RG |head -n 1 | tr "\t" " ")
if [[ -z "${RG}" ]]; then
	SM=$(echo ${TAG} | cut -f 2 -d "_")
	LB=$(echo ${TAG} | cut -f 4 -d "_")
	PL="ILLUMINA"
	PU="unit1"
	ID=${SM}_$(date +"%Y%m%d")_giraffe
	RG="${ID} ${LB} ${SM} ${PL} ${PU}"
fi

##vg mapping
#copy ref files
WORKDIR=$(dirname ${OUTPUTTAG})/work/
check_mkdir ${WORKDIR}
tmpref=${WORKDIR}/$(basename ${reference})
tmpdist=${WORKDIR}/$(basename ${dist})
tmpmin=${WORKDIR}/$(basename ${minimizer})
cp -f ${reference} ${WORKDIR}/ || exit $?
cp -f ${dist} ${WORKDIR}/ || exit $?
cp -f ${minimizer} ${WORKDIR}/ || exit $?
is_file_exists ${tmpref}
is_file_exists ${tmpdist}
is_file_exists ${tmpmin}
if [[ ${reference} =~ .gg$ ]]; then
	tmpgbwt=${WORKDIR}/$(basename ${gbwt})
	tmpxg=${WORKDIR}/$(basename ${xg})
	cp -f ${gbwt} ${WORKDIR}/ || exit $?
	cp -f ${xg} ${WORKDIR}/ || exit $?
	is_file_exists ${tmpgbwt}
	is_file_exists ${tmpxg}
fi
chmod -R 770 ${WORKDIR}

if [[ ${reference} =~ .gg$ ]]; then
	vg giraffe --progress --read-group "${RG}" --sample ${TAG} -g ${tmpref} -H ${tmpgbwt} -m ${tmpmin} -d ${tmpdist} -f ${fastq1} -f ${fastq2} -t ${NSLOTS} -o gam > ${OUTPUTTAG}.vg.giraffe.gam || exit $?
else
	vg giraffe --progress --read-group "${RG}" --sample ${TAG} -Z ${tmpref} -m ${tmpmin} -d ${tmpdist} -f ${fastq1} -f ${fastq2} -t ${NSLOTS} -o gam > ${OUTPUTTAG}.vg.giraffe.gam || exit $?
fi
if [[ "${mkbam}" == "TRUE" ]]; then
	if [[ ${reference} =~ .gg$ ]]; then
		pathref=${xg}
	else
		pathref=${tmpref}
	fi
	#extract path list
	vg paths -x ${pathref} -M | awk -F "\t" '{if($2=="REFERENCE"||$2=="GENERIC"){print $1}}' | awk -F '[\\[\\]]' '{print $1}' | sort | uniq > ${OUTPUTTAG}.surject.paths.txt
	#check reference
	junk=$(cat ${OUTPUTTAG}.surject.paths.txt | cut -f 1 -d "#" | sort | uniq | wc -l)
	junkreftag=$(cat ${OUTPUTTAG}.surject.paths.txt | grep "GRCh38")
	if [[ "${junk}" -gt 1 ]]; then
		if [[ -n "${junkreftag}" ]]; then
			mv ${OUTPUTTAG}.surject.paths.txt ${OUTPUTTAG}.surject.paths.txt.bk
			cat ${OUTPUTTAG}.surject.paths.txt.bk | grep "GRCh38" > ${OUTPUTTAG}.surject.paths.txt || exit $?
			rm ${OUTPUTTAG}.surject.paths.txt.bk
		else
			echo "Path file: ${OUTPUTTAG}.surject.paths.txt has more than one reference, but could not find GRCh38 path."
			exit 1
		fi
	fi
	vg surject -x ${pathref} --read-group "${RG}" --sample ${TAG} --interleaved --max-frag-len ${maxlen} --bam-output --prune-low-cplx ${OUTPUTTAG}.vg.giraffe.gam > ${OUTPUTTAG}.vg.giraffe.bam || exit $?
fi

##clean files
if [[ "${rmfastq}" == "TRUE" ]] ; then
	rm ${fastq1} ${fastq2}
fi
if [[ -f "${tmpref}" ]]; then
	rm ${tmpref}
fi
if [[ -f "${tmpdist}" ]]; then
	rm ${tmpdist}
fi
if [[ -f "${tmpmin}" ]]; then
	rm ${tmpmin}
fi
if [[ -f "${tmpgbwt}" ]]; then
	rm ${tmpgbwt}
fi
if [[ -f "${tmpxg}" ]]; then
	rm ${tmpxg}
fi
rmdir ${WORKDIR}
