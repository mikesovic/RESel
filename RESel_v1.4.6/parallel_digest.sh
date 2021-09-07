#PBS -l walltime=08:00:00
#PBS -l nodes=1:ppn=2
#PBS -A PAS1792
#PBS -j oe

#argument 1 is fasta genome file to operate on
#argument 2 is name of RE_info file (tab delimited)
#argument 3 is optional and is a bed-formatted (4 cols) features file.
#these arguments are passed in from the wrapper script.

cd $PBS_O_WORKDIR

name=`echo "${1}" | sed -E 's/.fasta$//' | sed -E 's/.fa$//' | sed -E 's/.fna$//'`

cut -f1 ${2} > seqs_${2}

REs1=`cut -f2 ${2} | sort -u | tr "\n" "_"`
REs=`echo ${REs1} | sed 's/_$//'`

i=`grep -c '^>' ${1}`
lengthi=`echo ${i} | wc -c`

mkdir ${name}_out/${REs}
touch ${name}_out/${REs}/cutsites.tsv

for j in $(seq 1 ${i})
do
	lengthj=`echo ${j} | wc -c`
	zerosNeeded="$(( ${lengthi} - ${lengthj} ))"
	
	if [[ ${zerosNeeded} == 1 ]] ; then 
	  j=`echo "0${j}"`
	fi
	
	if [[ ${zerosNeeded} == 2 ]] ; then 
	  j=`echo "00${j}"`
	fi
	
	if [[ ${zerosNeeded} == 3 ]] ; then 
	  j=`echo "000${j}"`
	fi
	
	if [[ ${zerosNeeded} == 4 ]] ; then 
	  j=`echo "0000${j}"`
	fi
	
	if [[ ${zerosNeeded} == 5 ]] ; then 
	  j=`echo "00000${j}"`
	fi
	
	if [[ ${zerosNeeded} == 6 ]] ; then 
	  j=`echo "000000${j}"`
	fi

chr=`grep '>' unwrapped/${name}_chr${j}_unwrapped.fa | sed 's/>//'`
grep -F -f seqs_${2} -aob uppers/${name}_${j}.txt | awk -v chr=${chr} -F ":" 'BEGIN {OFS="\t"} {print $1,$2,chr}' >> ${name}_out/${REs}/cutsites.tsv
done

mv ${2} ${name}_out/${REs}/${2}
rm seqs_${2}

module load R/4.0.2-gnu9.1

if [ ${3} ]
 then Rscript --verbose ~/software/RESelv1.4.6/scripts/analyze_dig.R ${name} ${REs} ${name}_out/${REs}/${2} ${3} > ${name}_out/${REs}/digest.Rout
 else  Rscript --verbose ~/software/RESelv1.4.6/scripts/analyze_dig.R ${name} ${REs} ${name}_out/${REs}/${2} > ${name}_out/${REs}/digest.Rout
fi
