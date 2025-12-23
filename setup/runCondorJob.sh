#!/bin/bash
source /cvmfs/alice.cern.ch/etc/login.sh
eval $(alienv printenv VO_ALICE@pythia::v8304-23,VO_ALICE@ROOT::v6-26-10-alice5-2)

NEVENTS=10000
NCPUS=1

export INPUT_FILES_DIR=/user/panosch/UM/Projects/Condor/Code
export LOG_FILES_DIR=/data/alice/pchrist/UM/Project3000

outdir=/dcache/alice/panosch/UM/Project3000/Monash/job_$1

mkdir -p ${outdir}
mkdir -p ${LOG_FILES_DIR}

export WORKINGDIR=/data/alice/pchrist/UM/Project3000/test/
cd ${WORKINGDIR}
echo ${WORKINGDIR}
cp ${INPUT_FILES_DIR}/* ${WORKINGDIR}

pwd
ls

##run the simulation                                                                                          
echo "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
source compile.sh
source generate.sh
#source analyse.sh
#aliroot -q -b fastCocktailProduction.C\($NEVENTS\)
#if [ ! -f galice.root ]; then
#        echo "Failed, trying again"
#	aliroot -q -b fastCocktailProduction.C\($NEVENTS\)
#fi
#if [ ! -f galice.root ]; then
#        echo "Failed yet again. Second time's the charm"
#	aliroot -q -b fastCocktailProduction.C\($NEVENTS\)
#fi
#if [ ! -f galice.root ]; then
#        echo "Failed yet again. Third time's the charm"
#	aliroot -q -b fastCocktailProduction.C\($NEVENTS\)
#fi
#if [ ! -f galice.root ]; then
#echo "Forget it. Clean tree output up, this is bogus, sorry"
#        #rm *.root
#fi
echo "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"

#mv ${LOG_FILES_DIR}/job$1.error ${outdir}/
#mv ${LOG_FILES_DIR}/job$1.log ${outdir}/
#mv ${LOG_FILES_DIR}/job$1.out ${outdir}/
mv ${WORKINGDIR}/*.root ${outdir}/

echo "Anything left?"
#rm -rf ${WORKINGDIR}
#rm -rf ${LOG_FILES_DIR}/job$1.error
#rm -rf ${LOG_FILES_DIR}/job$1.log
#rm -rf ${LOG_FILES_DIR}/job$1.out
ls
echo "Done! Enjoy!"
