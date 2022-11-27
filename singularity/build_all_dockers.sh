#!/bin/bash
TWD=`pwd` ;

echo -ne "Running script in ${TWD}\n\n\n" ; 
for DF in `find ${TWD} -name "Dockerfile" | grep -Pv 'CircTools' |  grep -Pv 'FUCH' `; do
	echo "Found Dockerfile ${DF}" ; 
	ID=`dirname ${DF}` ; 
	/bin/echo -ne "\tIts directory : ${ID}\n" ;
	/bin/echo -ne "\n\n\n\n" ; 
	ID_BN=`basename ${ID}`; 
	DOCKERNAME=`/bin/echo -ne "local/${ID_BN}" | tr '[:upper:]' '[:lower:]'`;
	echo "To use dockername ${DOCKERNAME}" ; 
	cd ${ID} && docker build . -t ${DOCKERNAME}  2>&1| tail --lines=3 && cd ${TWD} ; 
	/bin/echo -ne "\n\n\n\n" ; 
done ;



