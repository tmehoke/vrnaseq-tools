
#===================================================================================================
#
# This script will eventually be fully automated.
# For now it is a record of I have done manually, moving in the direction of automation.
# Each line should be sequentially understood and then manually run if it makes sense to do so.
#
#===================================================================================================


# create directory structure for kraken database
KRAKEN_HOME="/data/indices/kraken"
DB="fludb_20150325"
BASE="$KRAKEN_HOME/$DB"
mkdir -p $KRAKEN_HOME/$DB/{taxonomy,library,raw/individuals/bad}


# creating a kraken flu database using data downloaded from:
#   http://www.fludb.org/brc/influenza_sequence_search_segment_display.spg?method=ShowCleanSearch&decorator=influenza

# From this site, two searches were made, one for Influenza A, one for Influenza B

#	Data to Return: Segment / Nucleotide
#	Virus Type: A
#	Sub Type: (blank)
#	Strain Name: (blank)
#	Select Segments: All
#	Complete Sequences: Include Partial Sequences
#	2009 pH1N1 Sequences: Include pH1N1 sequences
#	Date Range: (blank)
#	Host: All
#	Geographic Grouping: All
#	Country: (blank)

# numbers of results
# 1) PB2 = 32,145
# 2) PB1 = 32,079
# 3) PA = 31,992
# 4) HA = 80,297
# 5) NP = 32,335
# 6) NA = 51,986
# 7) M = 43,110
# 8) NS = 33,421

# Total: 337,365

#	Data to Return: Segment / Nucleotide
#	Virus Type: B
#	Sub Type: (blank)
#	Strain Name: (blank)
#	Select Segments: All
#	Complete Sequences: Include Partial Sequences
#	2009 pH1N1 Sequences: Include pH1N1 sequences
#	Date Range: (blank)
#	Host: All
#	Geographic Grouping: All
#	Country: (blank)

# numbers of results
# 1) PB1 = 2,073
# 2) PB2 = 2,074
# 3) PA = 2,057
# 4) HA = 7,998
# 5) NP = 2,092
# 6) NA/NB = 5,594
# 7) M1/BM2 = 2,106
# 8) NS = 3,221

# Total: 27,215

# All data from each was saved as two FASTA files: InfA.fasta () and InfB.fasta ()
# The FASTA headers are currently formatted as follows:
# >gb:K00429|Organism:Influenza A virus A/seal/Mass/1/1980|Segment:4|Subtype:H7N7|Host:Sea Mammal


# I also want to have the human, pig, dog, chicken genomes, as well as the UniVec database, and potentially also all of RefSeq complete


# The goal is to have a taxonomy where:
#	the first level is the segment
#	the second level is the segment type (for HA and NA), or the overall serotype (for the other segments)
#	the third level is the year
#	the fourth level is the specific strain itself

# It would be great to figure out which of these values is the most important to account for flu diversity

# I think I just need to create several of these databases and see what they look like


# There is an issue with nine files where the host is labeled as "Host:	Tufted Duck" or "Host:	Northern Shoveler" with a tab in the middle of the header string
grep $'\t' "$BASE/raw/InfA.fasta"
grep $'\t' "$BASE/raw/InfB.fasta" # there are no tabs in InfB

# fix these by removing the tab (there should be no tabs initially present in this FASTA file)
sed -i 's/\t//' "$BASE/raw/InfA.fasta"
#sed -i 's/\t//' "$BASE/raw/InfB.fasta" # not needed

# remove the blank line between entries, removing line breaks in the nucleotide sequence
sed -i -e '/^$/d' -e ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' "$BASE/raw/InfA.fasta"
sed -i -e '/^$/d' -e ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' "$BASE/raw/InfB.fasta"

# remove any sequences that say they are from a 'mixed' subtype
sed '$!N;s/\n/\t/' "$BASE/raw/InfA.fasta" | grep -v "|Subtype:[mM]ixed|" | grep -v "|Subtype:H?" | sed 's/\t/\n/' > "$BASE/raw/InfA.fixed.fasta"
sed '$!N;s/\n/\t/' "$BASE/raw/InfB.fasta" | grep -v "|Subtype:[mM]ixed|" | grep -v "|Subtype:H?" | sed 's/\t/\n/' > "$BASE/raw/InfB.fixed.fasta"
# InfB only has two sequences that are mixed, and none that are H?

# list all hemagglutinin present
sed '$!N;s/\n/\t/' "$BASE/raw/InfA.fixed.fasta" | cut -f1 | cut -d"|" -f4 | cut -d":" -f2 | egrep '^H([0-9]*)N([0-9]*$)' | sed 's/H\([0-9]*\)N[0-9]*/\1/' | sort -n | uniq

# list all neuraminidase present
sed '$!N;s/\n/\t/' "$BASE/raw/InfA.fixed.fasta" | cut -f1 | cut -d"|" -f4 | cut -d":" -f2 | egrep '^H([0-9]*)N([0-9]*$)' | sed 's/H[0-9]*N\([0-9]*\)/\1/' | sort -n | uniq

# list all serotypes present
sed '$!N;s/\n/\t/' "$BASE/raw/InfA.fixed.fasta" | cut -f1 | cut -d"|" -f4 | cut -d":" -f2 | egrep '^H([0-9]*)N([0-9]*$)' | sort | uniq
sed '$!N;s/\n/\t/' "$BASE/raw/InfB.fixed.fasta" | cut -f1 | cut -d"|" -f4 | cut -d":" -f2 | egrep '^H([0-9]*)N([0-9]*$)' | sort | uniq
# Influenza B doesn't have a subtype...so I'm not sure what to do

# there are two samples that are listed as H4N48 for Influenza A virus A/mallard/Minnesota/33/2000
# this is actually an H4N4,8 (which I'm not sure what that means)

#---------------------------------------------------------------------------------------------------

# initializing names.dmp and nodes.dmp
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz -O "$KRAKEN_HOME/$DB/taxonomy/taxdump.tar.gz"
tar -xzf "$KRAKEN_HOME/$DB/taxonomy/taxdump.tar.gz"

#---------------------------------------------------------------------------------------------------

# initializing names.new and nodes.new

# new taxids will start at this value + 1
offset=1000000000

# initialize names.new with Influenza A segments
echo -e "$(($offset+1))\t|\tPB2\t|\t\t|\tscientific name\t|" >> "$BASE/taxonomy/names.new"
echo -e "$(($offset+2))\t|\tPB1\t|\t\t|\tscientific name\t|" >> "$BASE/taxonomy/names.new"
echo -e "$(($offset+3))\t|\tPA,PA-X\t|\t\t|\tscientific name\t|" >> "$BASE/taxonomy/names.new"
echo -e "$(($offset+4))\t|\tHA\t|\t\t|\tscientific name\t|" >> "$BASE/taxonomy/names.new"
echo -e "$(($offset+5))\t|\tNP\t|\t\t|\tscientific name\t|" >> "$BASE/taxonomy/names.new"
echo -e "$(($offset+6))\t|\tNA\t|\t\t|\tscientific name\t|" >> "$BASE/taxonomy/names.new"
echo -e "$(($offset+7))\t|\tM1,M2\t|\t\t|\tscientific name\t|" >> "$BASE/taxonomy/names.new"
echo -e "$(($offset+8))\t|\tNS1,NS2\t|\t\t|\tscientific name\t|" >> "$BASE/taxonomy/names.new"

# initialize names.new with Influenza B segments
echo -e "$(($offset+9))\t|\tPB1\t|\t\t|\tscientific name\t|" >> "$BASE/taxonomy/names.new"
echo -e "$(($offset+10))\t|\tPB2\t|\t\t|\tscientific name\t|" >> "$BASE/taxonomy/names.new"
echo -e "$(($offset+11))\t|\tPA\t|\t\t|\tscientific name\t|" >> "$BASE/taxonomy/names.new"
echo -e "$(($offset+12))\t|\tHA\t|\t\t|\tscientific name\t|" >> "$BASE/taxonomy/names.new"
echo -e "$(($offset+13))\t|\tNP\t|\t\t|\tscientific name\t|" >> "$BASE/taxonomy/names.new"
echo -e "$(($offset+14))\t|\tNA,NB\t|\t\t|\tscientific name\t|" >> "$BASE/taxonomy/names.new"
echo -e "$(($offset+15))\t|\tM1,BM2\t|\t\t|\tscientific name\t|" >> "$BASE/taxonomy/names.new"
echo -e "$(($offset+16))\t|\tNS\t|\t\t|\tscientific name\t|" >> "$BASE/taxonomy/names.new"

# initialize nodes.new with Influenza A segments
for i in {1..8}; do
	echo -e "$(($offset+$i))\t|\t11320\t|\tno rank\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|" >> "$BASE/taxonomy/nodes.new"
done
# initialize nodes.new with Influenza B segments
for i in {9..16}; do
	echo -e "$(($offset+$i))\t|\t11520\t|\tno rank\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|" >> "$BASE/taxonomy/nodes.new"
done

#---------------------------------------------------------------------------------------------------

# putting every year from 1918 to 2118 under every serotype present
for segment in {1..8}; do
	if [[ $segment -ne 4 && $segment -ne 6 ]]; then
		sed '$!N;s/\n/\t/' "$BASE/raw/InfA.fixed.fasta" | cut -f1 | grep "Segment:$segment" | \
				cut -d"|" -f4 | cut -d":" -f2 | egrep '^H([0-9]*)N([0-9]*$)' | sort | uniq | while read serotype; do
			taxid=$((1 + $(tail -n1 "$BASE/taxonomy/names.new" | cut -f1) ))
			echo -e "$taxid\t|\t$serotype\t|\t\t|\tscientific name\t|" >> "$BASE/taxonomy/names.new"
			echo -e "$taxid\t|\t$((segment+$offset))\t|\tno rank\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|" >> "$BASE/taxonomy/nodes.new"
			seq 1918 2118 | while read year; do
				newtaxid=$((1 + $(tail -n1 "$BASE/taxonomy/names.new" | cut -f1) ))
				echo -e "$newtaxid\t|\t$year\t|\t\t|\tscientific name\t|" >> "$BASE/taxonomy/names.new"
				echo -e "$newtaxid\t|\t$taxid\t|\tno rank\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|" >> "$BASE/taxonomy/nodes.new"
			done
		done
	elif [[ $segment -eq 4 ]]; then
		sed '$!N;s/\n/\t/' "$BASE/raw/InfA.fixed.fasta" | cut -f1 | grep "Segment:$segment" | \
				cut -d"|" -f4 | cut -d":" -f2 | egrep '^H([0-9]*)N([0-9]*$)' | sed 's/\(H[0-9]*\)N[0-9]*/\1/' | sort | uniq | while read serotype; do
			taxid=$((1 + $(tail -n1 "$BASE/taxonomy/names.new" | cut -f1) ))
			echo -e "$taxid\t|\t$serotype\t|\t\t|\tscientific name\t|" >> "$BASE/taxonomy/names.new"
			echo -e "$taxid\t|\t$((segment+$offset))\t|\tno rank\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|" >> "$BASE/taxonomy/nodes.new"
			seq 1918 2118 | while read year; do
				newtaxid=$((1 + $(tail -n1 "$BASE/taxonomy/names.new" | cut -f1) ))
				echo -e "$newtaxid\t|\t$year\t|\t\t|\tscientific name\t|" >> "$BASE/taxonomy/names.new"
				echo -e "$newtaxid\t|\t$taxid\t|\tno rank\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|" >> "$BASE/taxonomy/nodes.new"
			done
		done
	elif [[ $segment -eq 6 ]]; then
		sed '$!N;s/\n/\t/' "$BASE/raw/InfA.fixed.fasta" | cut -f1 | grep "Segment:$segment" | \
				cut -d"|" -f4 | cut -d":" -f2 | egrep '^H([0-9]*)N([0-9]*$)' | sed 's/H[0-9]*\(N[0-9]*\)/\1/' | sort | uniq | while read serotype; do
			taxid=$((1 + $(tail -n1 "$BASE/taxonomy/names.new" | cut -f1) ))
			echo -e "$taxid\t|\t$serotype\t|\t\t|\tscientific name\t|" >> "$BASE/taxonomy/names.new"
			echo -e "$taxid\t|\t$((segment+$offset))\t|\tno rank\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|" >> "$BASE/taxonomy/nodes.new"
			seq 1918 2118 | while read year; do
				newtaxid=$((1 + $(tail -n1 "$BASE/taxonomy/names.new" | cut -f1) ))
				echo -e "$newtaxid\t|\t$year\t|\t\t|\tscientific name\t|" >> "$BASE/taxonomy/names.new"
				echo -e "$newtaxid\t|\t$taxid\t|\tno rank\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|" >> "$BASE/taxonomy/nodes.new"
			done
		done
	fi
done

# create .root files that only contain the potential parent nodes - for future searching
cp "$BASE/taxonomy/names.new" "$BASE/taxonomy/names.root"
cp "$BASE/taxonomy/nodes.new" "$BASE/taxonomy/nodes.root"

#---------------------------------------------------------------------------------------------------

# format FASTA to be single files, and add GI numbers

# pull most recent livelist from NCBI
wget ftp://ftp.ncbi.nih.gov/genbank/livelists/GbAccList.0412.2015.gz -O "$KRAKEN_HOME/DB/GbAccList.0412.2015.gz"
gzip -d "$KRAKEN_HOME/DB/GbAccList.0412.2015.gz"

# pull GI list for all Influenza samples from NCBI
# How I actually created this "$gilist" file was to search "Influenza" on NCBI and then save all gi numbers as "sequence.gi.txt"
# I am going to figure out how to automatically perform this search from the script
utils="http://eutils.ncbi.nlm.nih.gov/entrez/eutils"
db="nuccore"
term="Influenza"
wget -q "$utils/esearch.fcgi?db=$db&term=$term&usehistory=y&tool=efetch" -O temp.txt


# FIX THIS LINE
gilist="/path/to/search/results"


# sort these two lists and then join them
sort -t "," -k3,3 "$KRAKEN_HOME/DB/GbAccList.0412.2015.gz" > "$BASE/raw/GbAccList.0412.2015.sorted"
sort $gilist > "$BASE/raw/sequence.gi.txt.sorted"
join -1 1 -2 3 -t "," "$BASE/raw/sequence.gi.txt.sorted" "$BASE/raw/GbAccList.0412.2015.sorted" > "$BASE/raw/gb2gi.txt"
paste -d"," \
	<(cut -d"," -f2-3 gb2gi.txt) \
	<(cut -d"," -f1 gb2gi.txt) > "$BASE/raw/gb2gi.txt.temp" \
	&& mv "$BASE/raw/gb2gi.txt.temp" "$BASE/raw/gb2gi.txt"

#=================================================
# The below commands were run as a Torque job

# add GI numbers to the beginning of each file
acclist="$BASE/raw/gb2gi.txt"
gb2gi() {
	header="$1"
	seq="$2"
	acc=$(echo $header | cut -d"|" -f2)
	if [[ ! -s "$BASE/raw/individuals/$acc.fna" ]]; then
		gi=$(grep "^$acc," "$acclist" | cut -d"|" -f2)
		if [[ -n $gi ]]; then
			echo $header | sed "s/>/>gi|$gi|/" > "$BASE/raw/individuals/$acc.fna"
			echo $seq >> "$BASE/raw/individuals/$acc.fna"
		else
			echo "Warning: $acc not found in NCBI Influenza A list, checking full GenBank Accession list" >&2
			gi=$(grep "^$acc," "/data/indices/ftp.ncbi.nih.gov/genbank/livelists/GbAccList.0412.2015" | cut -d"," -f3)
			if [[ -n $gi ]]; then
				echo $header | sed "s/>/>gi|$gi|/" > "$BASE/raw/individuals/$acc.fna"
				echo $seq >> "$BASE/raw/individuals/$acc.fna"
				else
				echo "Warning: $acc not found" >&2
			fi
		fi
	fi
}

# export variables for GNU parallel
export BASE
export acclist
export -f gb2gi

parallel -d $'\n' --colsep $'\t' -a <(sed '$!N;s/\n/\t/' "$BASE/raw/InfA.fixed.fasta") -n 1 -P 16 gb2gi {1} {2}

# endo of gb2gi torque job
#=================================================

#=================================================
# The below commands were run as a Torque job

taxonomy_build() {

	# parse input arguments
	newtaxid="$1"
	fn="$2"

	# pull out relevant information from header
	segment=$(head -n1 "$fn" | cut -d"|" -f6 | cut -d":" -f2)
	serotype=$(head -n1 "$fn" | cut -d"|" -f7 | cut -d":" -f2 | egrep '^H[0-9]+N[0-9]+$')
	year=$(head -n1 "$fn" | cut -d"|" -f5 | egrep '.*/([0-9]+)$' | sed 's/.*\/\([0-9]*\)/\1/' | sed 's/\(^0[0-9]$\)/20\1/' | \
		sed 's/\(^[2-9][0-9]$\)/19\1/' | sed 's/\(^1[0-4]$\)/20\1/' | sed 's/\(^1[5-9]$\)/19\1/')

	# remove files that don't contain good information
	if [[ -z "$segment" || -z "$serotype" || -z "$year" ]]; then
		mv "$fn" "$BASE/raw/individuals/bad/"
	else

		# create gi_taxid_nucl.new - two column file: gi, taxid
		head -n1 "$fn" | cut -d"|" -f2 | sed 's/[A-Z]*\([0-9]*\)/\1/' | nl -nrz | \
			awk -v TAXID="$newtaxid" '{
				printf "%d\t%d\n", $2, TAXID
			}' >> "$BASE/taxonomy/gi_taxid_nucl.new"

		# subset serotype if on an HA or NA segment
		if [[ "$segment" -eq 4 ]]; then
			serotype="${serotype%%N*}"
		elif [[ "$segment" -eq 6 ]]; then
			serotype="N${serotype#H[0-9]*N}"
		fi

		# add to names.new
		head -n1 "$fn" | \
			awk -F"|" -v TAXID="$newtaxid" -v SEGMENT="$segment" -v YEAR="$year" '{
			ORG=substr($5,10);
			TYPE=substr($7,9);
			printf "%d\t|\t%s (%s)\t|\t\t|\tscientific name\t|\n%d\t|\t%s\t|\t\t|\tsegment\t|\n%d\t|\t%d\t|\t\t|\tsubtype\t|\n%d\t|\t%d\t|\t\t|\tyear\t|\n", NR+TAXID-1, ORG, TYPE, NR+TAXID-1, SEGMENT, NR+TAXID-1, TYPE, NR+TAXID-1, YEAR
			}' >> "$BASE/taxonomy/names.new"

		# identify parent ID
		rootid=$(join -j1 \
			<(awk -F"\t" -v SEGMENTID=$((segment+1)) '{if($3 == SEGMENTID){print $1}}' "$BASE/taxonomy/nodes.new" | sort -k1b,1) \
			"$BASE/taxonomy/names.root" | \
			grep "| $serotype |" | cut -d" " -f1)

		# add to nodes.new
		join -j1 \
				<(awk -F"\t" -v ROOTID="$rootid" '{if($3 == ROOTID){print $1}}' "$BASE/taxonomy/nodes.new" | sort -k1b,1) \
				"$BASE/taxonomy/names.root" | \
			grep "| $year |" | cut -d" " -f1 | \
			awk -v TAXID="$newtaxid" '{
				printf "%d\t|\t%d\t|\tno rank\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\n", TAXID, $0
			}' >> "$BASE/taxonomy/nodes.new"
	fi
}

# create list of new taxids + filenames
offset=2000000000
find "$BASE/raw/individuals/" -maxdepth 1 -name "*.fna" -print0 | \
	awk -v OFFSET=$offset 'BEGIN {RS="\000"};  {printf("%d\t%s\n", NR+OFFSET, $0)}' > "$BASE/library/filenames.txt"

# export variables for GNU parallel
export BASE
export -f taxonomy_build

# add all files to taxonomy
parallel -d $'\n' --colsep $'\t' -a "$BASE/library/filenames.txt" -n 1 -P 16 taxonomy_build {1} {2}


# end of taxonomy-build torque job
#=================================================

# append new taxonomy files to full NCBI taxonomy files
cat "$BASE/taxonomy/names.new" >> "$BASE/taxonomy/names.dmp"
cat "$BASE/taxonomy/nodes.new" >> "$BASE/taxonomy/nodes.dmp"
cat "$BASE/taxonomy/gi_taxid_nucl.new" >> "$BASE/taxonomy/gi_taxid_nucl.dmp"

#---------------------------------------------------------------------------------------------------

# add to library
find "$BASE/raw/individuals" -maxdepth 1 -name "*.fna" -print0 | while read -d $'\0' fn; do
	kraken-build --add-to-library "$fn" --db "$BASE/"
done

#---------------------------------------------------------------------------------------------------

# build library
kraken-build --build --db "$BASE" --threads 64 | tee "$BASE/kraken-build.log" &

#---------------------------------------------------------------------------------------------------

