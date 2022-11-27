#!/bin/bash
TWD=`pwd` ;

echo -ne "Running script in ${TWD}\n\n\n" ; 
for DF in `find ${TWD} -name "build_sing.sh"  `; do
	echo "Found singularity script ${DF}" ; 
	ID=`dirname ${DF}` ; 
	/bin/echo -ne "\tIts directory : ${ID}\n" ;
	/bin/echo -ne "\n\n\n\n" ; 
	ID_BN=`basename ${ID}`; 
	cd ${ID_BN} ; 
	./build_sing.sh
	cd ${TWD} ; 
	/bin/echo -ne "\n\n\n\n" ; 
done ;



