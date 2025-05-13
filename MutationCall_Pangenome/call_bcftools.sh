#! /bin/sh

samtoolsDIR=/usr/local/package/samtools/1.19/bin/
export PATH=${samtoolsDIR};${PATH}
blatDIR=/home/package/blat/linux.x86_64/
export PATH=${blatDIR};${PATH}
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

#clone from call_bcftools.v2.sh for github

fasta=/home/reference/Human_Pangenome_Reference_Consortium/Graphs/GRCh38/AF_filtered/hprc-v1.1-mc-grch38.d9.gbz.GRCh38referencepaths.fa
##setting
score_difference=5
window_size=200
exclude_sam_flags=1024
EBcall_mapq=0
EBcall_baseq=15
base_qual_mpileup=15
Ncore=1

usage() {
	echo "Arguments=[tumor BAM] [tumor TAG] [normal bam] [(-d | --divconfig  [ARG])] [(-o | --output [ARG])] [(-c | --genomonconfig  [ARG])] [(-j | --jobfile [ARG])][(-J)]"
	echo "Options:"
	echo "  -h, --help"
	echo "  -f | --fasta [ARG] (option), fasta file for cram input, default = ${fasta}"
	echo "  -o | --output [ARG] (option), output path. default = ${DIR}/[tumor TAG]/[tumor TAG]"
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
		'--rm' )
			rminput="TRUE"
			shift 1
			;;
		'-o'|'--output' )
			if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
				echo "option requires an argument $1" 1>&2
				exit 1
			fi
			OUTPUTTAG=$2
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

##bcftools call
bcftools mpileup --threads ${Ncore} -f ${fasta} ${TBAM} ${NBAM} --annotate FORMAT/AD,FORMAT/DP,FORMAT/ADF,FORMAT/ADR,FORMAT/SP,INFO/AD,INFO/SCR --ff 4 -Q ${base_qual_mpileup} | bcftools call -f GQ --threads ${Ncore} -A -m -Oz -o ${OUTPUTTAG}.U1snRNA.bcfcall.vcf.gz || exit $?
tabix -f ${OUTPUTTAG}.U1snRNA.bcfcall.vcf.gz

bcftools norm -m-both ${OUTPUTTAG}.U1snRNA.bcfcall.vcf.gz| bcftools norm -f ${fasta} -Oz -o ${OUTPUTTAG}.BQ${base_qual_mpileup}.U1snRNA.bcfcall.lnorm.vcf.gz || exit $?
tabix -f ${OUTPUTTAG}.BQ${base_qual_mpileup}.U1snRNA.bcfcall.lnorm.vcf.gz
bcftools view ${OUTPUTTAG}.BQ${base_qual_mpileup}.U1snRNA.bcfcall.lnorm.vcf.gz -i 'FORMAT/AD[0:1] >= 3 && FORMAT/AD[1:1] <= 3' -Oz -o ${OUTPUTTAG}.U1snRNA.bcfcall.lnorm.3.vcf.gz || exit $?
rm ${OUTPUTTAG}.U1snRNA.bcfcall.vcf.gz ${OUTPUTTAG}.U1snRNA.bcfcall.vcf.gz.tbi

#/home/package/annovar/2019Oct24/convert2annovar.pl can be downloaded from ANNOVAR https://annovar.openbioinformatics.org/en/latest/#annovar-documentation
/home/package/annovar/2019Oct24/convert2annovar.pl -format vcf4old --allallele ${OUTPUTTAG}.U1snRNA.bcfcall.lnorm.3.vcf.gz > ${OUTPUTTAG}.U1snRNA.bcfcall.lnorm.3.tsv || exit $?
rm ${OUTPUTTAG}.U1snRNA.bcfcall.lnorm.3.vcf.gz

mutfilter realignment -s ${score_difference} -w ${window_size} -t ${OUTPUTTAG}.U1snRNA.bcfcall.lnorm.3.tsv -1 ${TBAM} -2 ${NBAM} -o ${OUTPUTTAG}.U1snRNA.bcfcall.lnorm.3.realign.tsv -r ${fasta} -b $(which blat) --exclude_sam_flags ${exclude_sam_flags} || exit $?
rm ${OUTPUTTAG}.U1snRNA.bcfcall.lnorm.3.tsv

EBFilter -f anno -q ${EBcall_mapq} -Q ${EBcall_baseq} --ff 4 ${OUTPUTTAG}.U1snRNA.bcfcall.lnorm.3.realign.tsv ${TBAM} ${Panel} ${OUTPUTTAG}.BQ${base_qual_mpileup}.U1snRNA.bcfcall.lnorm.3.EB.tsv --loption || exit $?
rm ${OUTPUTTAG}.U1snRNA.bcfcall.lnorm.3.realign.tsv
