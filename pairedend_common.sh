#!/bin/bash

# automatically exit the script on any error
set -e

runtime=$(date +"%Y%m%d%H%M%S%N")

# This script takes two paired end FASTQ files and pulls out reads common to both files
# such as what occurs after both are filtered separately

# 		by: Tom Mehoke
# 	written: September 25, 2013

# input: 2 FASTQ files
# output: 2 FASTQ files, each only containing reads common to both files
# output (optional): 2 FASTQ files containing the orphaned reads from each input

# requirements: none

usage()
{
cat << EOF
usage: $0 [options] -1 <forward.fastq> -2 <reverse.fastq>

This script filters a pair of paired-end FASTQ reads to only reads common to both files

OPTIONS:
   -h      show this message
   -u      report orphan reads to separate file (optional)
   -o      output path
   -s      output file suffix (default: .common)
   -r      orphan read suffix (default: .orphan)
   -1      forward reads FASTQ file
   -2      reverse reads FASTQ file
   -w      workdir (default: /tmp)

EOF
}

# set default values here
outpath="./"
suffix=".common"
orphan=".orphan"
unmatched=false
workdir="/tmp"

# parse input arguments
while getopts "huo:s:1:2:w:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
		u) unmatched=true ;;
		o) outpath=$OPTARG ;;
		s) suffix=$OPTARG ;;
		1) f=$OPTARG ;;
		2) r=$OPTARG ;;
		w) workdir=$OPTARG ;;
		?) usage; exit ;;
	esac
done

# pull off initial periods from suffix
suffix=${suffix/\./}
orphan=${orphan/\./}

# if necessary arguments are not present, display usage info and exit
if [[ -z "$f" ]] || [[ -z "$r" ]]; then
	 usage
	 exit 1
fi

#===================================================================================================

fbase=$(basename "$f")
rbase=$(basename "$r")

# converts FASTQ to 4-column tab-delimited file with one row per record, push spaced bit at end of the header to a new field
# use sort -c -k 1b,1 to first check if input is sorted - not sure what the actual test should be
sed '$!N;s/\n/\t/' "$f" | sed '$!N;s/\n/\t/' | sort -T "$workdir" -k 1b,1 | sed 's/ \([1-2]:N:0:[0-9]*\)/\t\1/g' > "$workdir/$fbase.format-$runtime"
sed '$!N;s/\n/\t/' "$r" | sed '$!N;s/\n/\t/' | sort -T "$workdir" -k 1b,1 | sed 's/ \([1-2]:N:0:[0-9]*\)/\t\1/g' > "$workdir/$rbase.format-$runtime"

# pull out matched reads - format back to normal paired-end FASTQ file
join -t $'\t' -1 1 -2 1 "$workdir/$fbase.format-$runtime" "$workdir/$rbase.format-$runtime" | \
		awk -F"\t" -v F="$fbase-$runtime" -v R="$rbase-$runtime" -v W="$workdir/" '{print $1" "$2"\n"$3"\n"$4"\n"$5 > W F ".common"; print $1" "$6"\n"$7"\n"$8"\n"$9 > W R ".common"}'
mv "$workdir/$fbase-$runtime.common" "$outpath/${fbase%.fastq}.$suffix.fastq"
mv "$workdir/$rbase-$runtime.common" "$outpath/${rbase%.fastq}.$suffix.fastq"

# pull out orphan reads if -u flag is present
if $unmatched; then
	join -t $'\t' -v 1 -1 1 -2 1 "$workdir/$fbase.format-$runtime" "$workdir/$rbase.format-$runtime" | \
			awk -F"\t" -v F="$fbase-$runtime" -v W="$workdir/" '{print $1" "$2"\n"$3"\n"$4"\n"$5 > W F ".orphan"}'
	join -t $'\t' -v 2 -1 1 -2 1 "$workdir/$fbase.format-$runtime" "$workdir/$rbase.format-$runtime" | \
			awk -F"\t" -v R="$rbase-$runtime" -v W="$workdir/" '{print $1" "$2"\n"$3"\n"$4"\n"$5 > W R ".orphan"}'
	mv "$workdir/$fbase-$runtime.orphan" "$outpath/$fbase%.fastq}.$orphan.fastq"
	mv "$workdir/$rbase-$runtime.orphan" "$outpath/$rbase%.fastq}.$orphan.fastq"
fi

# remove temporary files
rm "$workdir/$fbase.format-$runtime"
rm "$workdir/$rbase.format-$runtime"

#--eof--#

