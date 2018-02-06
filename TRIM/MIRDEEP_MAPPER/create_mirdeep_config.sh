#! /bin/sh

# Change the paths as necessary
# DATA_DIR is the location of raw fastq files
# PROJECT_DIR is the main directory for this analysis
# ALL_FILE_LIST is the list of renamed sample fastq files

export PROJECT_DIR="/home/smalkaram/Projects/travis20180117151534"
export DATA_DIR="/home/smalkaram/Projects/travis20180117151534/raw_data/file_links"
export ALL_FILE_LIST="/home/smalkaram/Projects/travis20180117151534/ALL_FILE_LIST"

# This script generates a mapping file with original sample names mapped to sample names compatible with miRDeep2 program
# Run as 7.create_mirdeep_config.sh ALL_FILE_LIST > mapping_file2.txt
# *** For an inidentified reason, the following names ending in ??x are not compatible with mirDeep2 program. So change them as indicated below
#changed aax to aza
#changed abx to azb

INPUT=$1

LIMIT=`cat $INPUT | wc -l`
NUM=1
CODE=""
for ONE in {a..z}
do
	for TWO in {a..z}
	do
		for TRI in {a..z}
		do
			if [ $NUM -gt $LIMIT ]; then
				break
			fi
			code=$ONE""$TWO""$TRI
			CODE="$CODE $code"
			NUM=`expr $NUM + 1`
		done
	done
done

CODE=($CODE)

NUM=0
exec<"$INPUT"
while read line
do
	L=($line)
	echo "${L[0]} ${CODE[$NUM]}"
	NUM=`expr $NUM + 1`
done
