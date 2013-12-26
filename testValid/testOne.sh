#!/bin/bash

# Use "-e" to exit if any command has rc != 0.
# This is critical to the functioning of this script.
set -e
##set -v

usage() {
  echo ''
  echo testOne.sh: Error: $*
  echo 'parms:  cifFile  refDir  workDir'
  echo ''
  echo 'CAUTION: The workDir will be deleted and recreated.'
  echo 'The workDir must start with "/" and must end with "temp.test"'
  echo ''
  echo 'The cifFile must start with "/" and end with ".cif"'
  echo ''
  echo 'The refDir must start with "/"'
  echo ''
  echo 'Example: ./testOne.sh ~/pylada/testValid/icsd_633029/icsd_633029.cif ~/peregrine/testlada.FeO.icsd_633029 /tmp/temp.test'
  echo ''
  echo 'Error: see above'
  exit 1
}

if [ $# -ne 3 ]; then
  usage Wrong num parms
fi

cifFile=$1
refDir=$2
workDir=$3

if [[ ! "$cifFile" =~ ^/.*\.cif$ ]]; then
  usage invalid cifFile $cifFile
fi
if [[ ! -e "$cifFile" ]]; then
  usage cifFile $cifFile not found
fi

if [[ ! "$refDir" =~ ^/ ]]; then
  usage invalid refDir $refDir
fi
if [[ ! -d "$refDir" ]]; then
  usage refDir $refDir not found
fi

if [[ ! "$workDir" =~ ^/.*/temp.test$ ]]; then
  usage invalid workDir $workDir
fi

rm -rf $workDir
mkdir $workDir
cd $workDir

cat ~/pylada/testValid/pylada.execMimic.json \
  | sed -e "s#XrefRootX#$refDir#g" -e "s#XworkRootX#$workDir#g" \
  > $workDir/pylada.execMimic.json

cp ~/pylada/test/highthroughput/{inputCif.py,test.py} .

mkdir structs
cp $cifFile structs

# Set bugLev 5
sed -i 's/setBugLev(0)/setBugLev(5)/' test.py

cat inputCif.py | sed -e "/^setTestValidProgram/ s|None|'~/nrelmat/nrelmat/execMimic.py'|" -e 's#/nopt/nrel/ecom/cid/vasp.pseudopotentials.a#/home/ssulliva/vladan/td.pseudos#' > inputTestValid.py

# diff inputCif.py inputTestValid.py

echo testDir: finished setup

ipython --TerminalIPythonApp.file_to_run='' << eofeofa > log.ipython.nonmag
import test
test.nonmagnetic_wave('pickle', inputpath='inputTestValid.py')
launch scattered --ppn 24 --account ss04 --queue batch --walltime=0:30:00
eofeofa

echo testDir: cwd: $(pwd)
mv log log.execMimic.nonmag
echo testDir: finished nonmag run

ipython --TerminalIPythonApp.file_to_run='' << eofeofb > log.ipython.mag
import test
explore pickle
test.magnetic_wave('pickle', inputpath='inputTestValid.py')
launch scattered --ppn 24 --account ss04 --queue short --walltime=0:30:00
eofeofb

mv log log.execMimic.mag
echo testDir: finished mag run

# All pbserr and stderr files must have 0 lines.
wc -c $(find . -name pbserr | sort)
wc -c $(find . -name stderr | sort)
wc -c $(find . -name pbserr | sort) | egrep '^0 total$' > /dev/null
wc -c $(find . -name stderr | sort) | egrep '^0 total$' > /dev/null
echo testDir: finished check of pbserr, stderr

# The execMimic logs must have no exc* or err* messages.
egrep -v 'exc|err' log.execMimic.nonmag > /dev/null
egrep -v 'exc|err' log.execMimic.mag    > /dev/null
echo testDir: finished check of execMimic logs


echo ''
echo testDir: all done


