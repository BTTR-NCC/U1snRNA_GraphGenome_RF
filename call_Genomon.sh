#! /bin/sh
#$ -S /bin/bash
#$ -e /home/ha6434/command/vg/log/
#$ -o /home/ha6434/command/vg/log/
#$ -cwd
#$ -l s_vmem=20G,mem_req=20G

source /home/ha6434/util/utility.sh
source /home/ha6434/util/source_share.sh

module use /usr/local/package/modulefiles/
export MODULEPATH=/home/ha6434/modulefiles:$MODULEPATH
module load samtools/1.9 blat/linux.x86_64

fasta=/home/ha6434/reference/Human_Pangenome_Reference_Consortium/Graphs/GRCh38/AF_filtered/hprc-v1.1-mc-grch38.d9.gbz.GRCh38referencepaths.fa
genomonconfig=/home/ha6434/command/multimapMut/config/Genomon2_multimap_somatic_setting.v2.config ;is_file_exists ${genomonconfig}

tsv2vcf=/home/ha6434/command/multimapMut/subscript/multimapcf.pair.php
bgzip=/home/ha6434/command/samtools/bgzip.sh
lfnorm=/home/ha6434/command/samtools/bcftools.leftnorm.sh
runR=/home/ha6434/command/util/run_R_3.5.0.sh
fisherfilt=/home/ha6434/command/multimapMut/subscript/somatic_fisher_filter.v2.R

source /home/ha6434/virturalenv3.7.17_OS8/Genomon2/bin/activate

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
source ${genomonconfig}

##Fisher call
fisher comparison -1 ${TBAM} -2 ${NBAM} -o ${OUTPUTTAG}.U1snRNA.Genomon.fisher.txt --ref_fa ${fasta} --samtools_path $(which samtools) --base_quality ${base_quality} --min_allele_freq ${min_allele_freq} --max_allele_freq ${max_allele_freq} --min_depth ${min_depth} --min_variant_read ${min_variant_read} --fisher_value ${fisher_value} --samtools_params "${samtools_params}" || exit $?

##Local realignment by blat
mutfilter realignment -s ${score_difference} -w ${window_size} -t ${OUTPUTTAG}.U1snRNA.Genomon.fisher.txt -1 ${TBAM} -2 ${NBAM} -o ${OUTPUTTAG}.U1snRNA.Genomon.realign.txt -r ${fasta} -b $(which blat) --exclude_sam_flags ${exclude_sam_flags} || exit $?
rm ${OUTPUTTAG}.U1snRNA.Genomon.fisher.txt|| exit $?

#Empirical Baysian mutation calling
deactivate
source /home/ha6434/virturalenv3.7.17_OS8/GenomonEBFilt/bin/activate

EBFilter -f anno -q ${EBcall_mapq} -Q ${EBcall_baseq} --ff 4 ${OUTPUTTAG}.U1snRNA.Genomon.realign.txt ${TBAM} ${Panel} ${OUTPUTTAG}.U1snRNA.Genomon.EB.tsv --loption || exit $?
rm ${OUTPUTTAG}.U1snRNA.Genomon.realign.txt || exit $?
