#!/bin/bash

# Added safety. If any command fails or a pipe breaks, the script will stop running.
set -eu -o pipefail -o verbose

# Usage: 
# bash ATCGmap_GenotypeCorrector.sh MyMethylome.ATCGmap.gz > MyMethylome.CGmap.corrected.gz

# Output format the same as CGmap (typical output from BSseeker2 or CGmapTools, but positions
# that show evidence for CpGs in the reads annotated as CH in the reference genome are
# annotated as CGnew. 

zcat $1 |awk 'BEGIN{OFS="\t";modifiedpos=0}{for (i=1; i<=NR; i++) {NR=i; line=$0; base=$2; pos1=$3; context=$4; \
	prevWatsonC=$8; prevWatsonT=$7; prevCrickG=$13; \
	prevWatsonNON_C=($6+$7+$9+$10); prevCrickNON_G=($11+$12+$14+$15); 			\
	getline; nextline=$0; pos2=$3; nextbase=$2; nextcontext=$4;\
	nextWatsonG=$9; nextCrickC=$14; nextCrickT=$11; \
	nextWatsonNON_G=($6+$7+$8+$10); nextCrickNON_C=($12+$13+$11+$15); 			\
	if  ( modifiedpos == pos1 ) {print oldline}	\
	else if ( base ~ /C/ && context ~ /CH/ && pos2 == ( pos1 + 1 ) && 				\
	( prevWatsonC >= 1 ) ) 			\
		{ if (	 \
			( nextWatsonNON_G + nextWatsonG ) > 5 &&  nextWatsonG >= nextWatsonNON_G )	\
	 		{sub(/CH[GH]/,"CGnew",line); sub(/--|CH[HG]|CG/,"CGnew",nextline);sub(/	(A|C|T)	/,"	G	",nextline); print line; oldline=nextline; modifiedpos=pos2}	\
	 	else if (  \
	 		( nextWatsonNON_G + nextWatsonG ) > 0 && ( nextCrickNON_C + nextCrickC ) > 0 && \
	 		(nextWatsonG/(nextWatsonNON_G + nextWatsonG )) > 0.2 && (nextCrickC/( nextCrickNON_C + nextCrickC )) >= 0.2 ) \
	 		{sub(/CH[GH]/,"CGnew",line) ; sub(/--|CH[HG]|CG/,"CGnew",nextline);sub(/	(A|C|T)	/,"	G	",nextline); print line; oldline=nextline; modifiedpos=pos2}	\
	 	else if (  \
	 		( nextWatsonNON_G + nextWatsonG ) == 0 && ( nextCrickNON_C + nextCrickC ) > 0 &&  	\
	 		nextCrickC/( nextCrickNON_C + nextCrickC ) >= 0.25 )								\
	 		{sub(/CH[GH]/,"CGnew",line) ; sub(/--|CH[HG]|CG/,"CGnew",nextline);sub(/	(A|C|T)	/,"	G	",nextline); print line; oldline=nextline; modifiedpos=pos2} \
	 	else if (  \
	 		( nextWatsonNON_G + nextWatsonG ) > 0 && ( nextCrickNON_C + nextCrickC ) == 0 &&  	\
	 		nextWatsonG >= nextWatsonNON_G )							\
	 		{sub(/CH[GH]/,"CGnew",line) ; sub(/--|CH[HG]|CG/,"CGnew",nextline);sub(/	(A|C|T)	/,"	G	",nextline); print line; oldline=nextline; modifiedpos=pos2} \
	 	else if (  \
	 		nextWatsonG > 1 && nextWatsonG > 0.2*nextWatsonNON_G && prevWatsonC/nextWatsonG >= 0.5 &&  	\
	 		prevWatsonC/nextWatsonG <= 1.5 )							\
	 		{sub(/CH[GH]/,"CGnew",line) ; sub(/--|CH[HG]|CG/,"CGnew",nextline);sub(/	(A|C|T)	/,"	G	",nextline); print line; oldline=nextline; modifiedpos=pos2} \
		else if (  \
	 		nextCrickC > 1 && prevCrickG > 0 && nextCrickC > 0.2*nextCrickNON_C && nextCrickC/prevCrickG >= 0.5 &&  	\
	 		nextCrickC/prevCrickG <= 1.5 )							\
	 		{sub(/CH[GH]/,"CGnew",line) ; sub(/--|CH[HG]|CG/,"CGnew",nextline);sub(/	(A|C|T)	/,"	G	",nextline); print line; oldline=nextline; modifiedpos=pos2} \
	  	else {print line}\
	 	} \
	else if (  \
	nextbase ~ /G/ && nextcontext ~ /CH/ && pos2 == ( pos1 + 1 ) && 			\
	( nextCrickC  >= 1 ) ) 		\
	 	{ if (	 \
			( prevCrickNON_G + prevCrickG ) > 5 && prevCrickG >= prevCrickNON_G ) \
	 		{sub(/--|CH[HG]|CG/,"CGnew",line);sub(/	(A|G|T)	/,"	C	",line); sub(/CH[GH]/,"CGnew",nextline) ; print line; oldline=nextline; modifiedpos=pos2}\
	 	else if (  \
	 		( prevCrickNON_G + prevCrickG ) > 0 && ( prevWatsonNON_C + prevWatsonC ) > 0 &&  \
	 		prevCrickG/( prevCrickNON_G + prevCrickG ) > 0.2 && prevWatsonC/( prevWatsonNON_C + prevWatsonC ) >= 0.2 ) \
	 		{sub(/--|CH[HG]|CG/,"CGnew",line);sub(/	(A|G|T)	/,"	C	",line); sub(/CH[GH]/,"CGnew",nextline) ; print line; oldline=nextline; modifiedpos=pos2} \
	 	else if (  \
	 		( prevCrickNON_G + prevCrickG ) == 0 && ( prevWatsonNON_C + prevWatsonC ) > 0 && 	\
	 		prevWatsonC/( prevWatsonNON_C + prevWatsonC ) >= 0.25 )								\
	 		{sub(/--|CH[HG]|CG/,"CGnew",line);sub(/	(A|G|T)	/,"	C	",line); sub(/CH[GH]/,"CGnew",nextline) ; print line; oldline=nextline; modifiedpos=pos2} \
	 	else if (  \
	 		( prevCrickNON_G + prevCrickG ) > 0 && ( prevWatsonNON_C + prevWatsonC ) == 0 &&  \
	 		 prevCrickG >= prevCrickNON_G ) \
	 		{sub(/--|CH[HG]|CG/,"CGnew",line);sub(/	(A|G|T)	/,"	C	",line); sub(/CH[GH]/,"CGnew",nextline) ; print line; oldline=nextline; modifiedpos=pos2} \
		else if (  \
	 		prevCrickG > 1 && prevCrickG > 0.2*prevCrickNON_G && nextCrickC/prevCrickG >= 0.5 &&  \
	 		 nextCrickC/prevCrickG <= 1.5 ) \
	 		{sub(/--|CH[HG]|CG/,"CGnew",line);sub(/	(A|G|T)	/,"	C	",line); sub(/CH[GH]/,"CGnew",nextline) ; print line; oldline=nextline; modifiedpos=pos2} \
		else if (  \
	 		prevWatsonC > 1 && nextWatsonG > 0 && prevWatsonC > 0.2*prevWatsonNON_C && prevWatsonC/nextWatsonG >= 0.5 &&  \
	 		 prevWatsonC/nextWatsonG <= 1.5 ) \
	 		{sub(/--|CH[HG]|CG/,"CGnew",line);sub(/	(A|G|T)	/,"	C	",line); sub(/CH[GH]/,"CGnew",nextline) ; print line; oldline=nextline; modifiedpos=pos2} \
	  	else {print line}\
	 	} \
	else{print line} }							\
	}' |awk 'BEGIN{OFS="\t"}{if ( $2 == "C" && ($8+$7) > 0 ){mC=$8;CT=($8+$7); print $1,$2,$3,$4,$5,mC/CT,mC,CT}else if( $2 == "G" && ($14+$11) > 0 ){mC=$14;CT=($14+$11); print $1,$2,$3,$4,$5,mC/CT,mC,CT}}'


