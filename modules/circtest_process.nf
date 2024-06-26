process circtest {

input:
    //cohort name
    val cohort_name

    //data from DCC
    path '*'

    //how to compare given cohort name
    path 'cohort_comp_conf'
    
    //plain or not
    val mode

output:
	file '*.circtest_results.zip'

shell:
'''
ls -alht
ls -alht /usr/local/bin/
##############################
#edit the files to have shorter names in the column headers; use the passed-in JSON file
for DCCF in LinearCount CircRNACount ; do
	python3 /usr/local/bin/name_edit.py !{cohort_comp_conf}  ${DCCF} s.${DCCF} ;
done ;
##############################
#distribute the columns into desired cohorts using
# an inline python3 script
COHORT_TO_RUN="!{cohort_name}" ; 
echo "COHORT_TO_RUN IS ${COHORT_TO_RUN}" ; 
COHORT_DIST="""
import os
import csv
import copy
import shutil
import json
linear_file='s.LinearCount'
circ_file='s.CircRNACount'
coord_file='CircCoordinates'

# read cohort/comparison data from JSON file
cohort_comp_conf_json_file='!{cohort_comp_conf}'
jreader=open(cohort_comp_conf_json_file)
cohort_comp_conf=json.load(jreader)
jreader.close()
cohort_comps=cohort_comp_conf['cohort_comps']
comp_names=cohort_comp_conf['comp_names']


def writeToFile(s,f):
    writer=open(f,'w')
    writer.write(s.strip()+'\\n')
    writer.close()
for comp in cohort_comps:
    if(not(comp=='${COHORT_TO_RUN}')):
        continue
    print(f'Processing for comparision {comp}')
    common_keys=['Chr','Start','End']
    os.makedirs(comp,exist_ok=True)
    shutil.copyfile(coord_file,comp+'/'+coord_file)
    tsv_files=[linear_file,circ_file]
    out_files=[comp+'/'+t for t in tsv_files]
    writeToFile(comp_names[comp][0],comp+'/G1.txt')
    writeToFile(comp_names[comp][1],comp+'/G2.txt')
    for file_pair in zip(tsv_files,out_files):
        with open(file_pair[0], mode='r') as csv_file:
            csv_reader = csv.DictReader(csv_file,delimiter='\\t')
            wroteHeader=False
            writer=open(file_pair[1],'w')
            recs_written=0
            for row in csv_reader:
                    if(not(wroteHeader)):
                        header=copy.deepcopy(common_keys)
                        header.extend(cohort_comps[comp])
                        wroteHeader=True
                        writer.write('\\t'.join(header)+'\\n')
                    row_keys_to_write=copy.deepcopy(common_keys)
                    row_keys_to_write.extend(cohort_comps[comp])
                    vals_to_write=[row[rk] for rk in row_keys_to_write]
                    writer.write('\\t'.join(vals_to_write)+'\\n')
                    recs_written=recs_written+1
            writer.close()

    #now write the counts to  analysis_col.txt
    counts_dict=cohort_comp_conf['comparison_counts']
    the_counts_this_comp=counts_dict[comp]
    val_to_write=the_counts_this_comp+3
    writer=open(comp+'/analysis_col.txt','w')
    writer.write(str(val_to_write))
    writer.close()    

"""
python3 -c "${COHORT_DIST}" ; 


BASE_DIR="${PWD}" ; 
for CT_DIR in `find */s.*|grep -Pi 'inear'|xargs -I {} dirname {}`; do
	cd ${BASE_DIR}/${CT_DIR} ; 
	DATA_IDX=`cat analysis_col.txt`; 
	NAME_A=`cat G1.txt`;
	NAME_B=`cat G2.txt`;
	#start a job
	set -x
    if [ "!{mode}" == "plain" ] ; then
        #plain uses only default parameters
	    circtest.cli.PLAIN.R s.CircRNACount s.LinearCount  CircCoordinates ${DATA_IDX} ${NAME_A} ${NAME_B}
    else 
        # this version of the script iterates over parameter combinations
        circtest.cli.R       s.CircRNACount s.LinearCount  CircCoordinates ${DATA_IDX} ${NAME_A} ${NAME_B}
    fi ;
	set +x
done ;
#zip up the cohort results
cd ${BASE_DIR} ;
DG1=`find */G1.txt|xargs dirname`; 
zip -r "${COHORT_TO_RUN}.circtest_results.zip" "${DG1}" ; 
'''

}


process circtest_plotting {


input:
	//zip file from circtest 
	path '*'
	//circatlas bed file for data joining
	path 'circatlas_bed.txt'
	//cohort comparison configuration
	path 'cohort_comp_conf'
	

output:
	path '*.plotting_results.zip'

shell:
'''
echo "Started at";
date
#unzip the file
unzip `find *.zip` ;
#EXTRACT COHORTNAME FROM ZIP
COHORT_NAME=`find *.zip | grep -Po '^[^\\.]+\\.'|tr -d "."`;
echo "From zip file, got cohort name ${COHORT_NAME}" ; 
#set up the output area
PLOTTING_DIR="${COHORT_NAME}_plottting" ;
mkdir -v ${PLOTTING_DIR}
#define the circatlas bed file
CA_BED="circatlas_bed.txt" ;
#run plotting
python3 /usr/local/bin/ct_aug.py !{params.plot_cutoff ? "-PVC "+params.plot_cutoff : ""} \
	 ${COHORT_NAME}	\
	 ${CA_BED}	\
	 ${PLOTTING_DIR}/${COHORT_NAME}.merged.tsv	\
	 ${PLOTTING_DIR} !{cohort_comp_conf} ; 
#zip up results of plotting
zip -r ${COHORT_NAME}.plotting_results.zip ${PLOTTING_DIR}
echo "Finished at";
date
'''

}


