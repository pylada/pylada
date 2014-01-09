#!/bin/bash

# Use "-e" to exit if any command has rc != 0.
# This is critical to the functioning of this script.
##xxxxxxxxxxset -e
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
  echo 'Example:'
  echo 'cd ~/pylada/testValid'
  echo './testOne.sh ~/pylada/testValid/icsd_633029/icsd_633029.cif ~/peregrine/testlada.FeO.icsd_633029 /tmp/temp.test'
  echo ''
  echo 'Error: see above'
  exit 1
}

# Quit if $1 is not 0.
checkRc() {
  rc=$1
  msg=$2
  if [[ "$rc" -ne 0 ]]; then
    echo Error: $msg
    exit 1
  fi
}

if [ $# -ne 3 ]; then
  usage Wrong num parms
fi

cifFile=$1
refDir=$2
workDir=$3

# Insure parameters exist
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

# Copy execMimic to workDir, inserting text for refDir, workDir.
cat ~/pylada/testValid/pylada.execMimic.json \
  | sed -e "s#XrefRootX#$refDir#g" -e "s#XworkRootX#$workDir#g" \
  > $workDir/pylada.execMimic.json
checkRc $? 'cannot setup execMimic'

# Set up the relaxation run
cp ~/pylada/test/highthroughput/inputCif.py .
checkRc $? 'cannot get inputCif.py'
cp ~/pylada/test/highthroughput/test.py .
checkRc $? 'cannot get test.py'

mkdir structs
cp $cifFile structs
checkRc $? 'cannot setup cif'

# Set bugLev 5
sed -i 's/setBugLev(0)/setBugLev(5)/' test.py
checkRc $? 'cannot mod test.py'

# Copy inputCif.py to inputTestR.py, setting the setTestValidProgram.
cat inputCif.py | sed -e "/^setTestValidProgram/ s|None|'~/nrelmat/nrelmat/execMimic.py'|" -e 's#/nopt/nrel/ecom/cid/vasp.pseudopotentials.a#/home/ssulliva/vladan/td.pseudos#' > inputTestR.py
checkRc $? 'cannot mod inputCif.py'

echo testOne: finished setup

echo ''
echo 'testOne: start non-magnetic relaxation'
echo ''

# Run non-magnetic tests
ipython --no-confirm-exit --colors=NoColor --TerminalIPythonApp.file_to_run='' << eofeofa > log.ipython.nonmag
import test
test.nonmagnetic_wave('pickle', inputpath='inputTestR.py')
launch scattered --ppn 24 --account ss04 --queue batch --walltime=0:30:00
eofeofa
checkRc $? 'ipython nonmag failed'

# Save the execMimic log
mv log log.execMimic.nonmag
checkRc $? 'cannot find nonmag log'
echo testOne: finished nonmag run





echo ''
echo 'testOne: start non-magnetic GW'
echo ''








# Set up the GW run
cp ~/pylada/test/highthroughput/inputGW.py .
checkRc $? 'cannot get inputGW.py'

# Copy inputGW.py to inputTestG.py, setting the setTestValidProgram.
cat inputGW.py | sed -e "/^setTestValidProgram/ s|None|'~/nrelmat/nrelmat/execMimic.py'|" -e 's#/nopt/nrel/ecom/cid/vasp.pseudopotentials.a#/home/ssulliva/vladan/td.pseudos#' > inputTestG.py
checkRc $? 'cannot mod inputGW.py'

echo testOne: finished setup

# Run the nonmag GW tests
ipython --no-confirm-exit --colors=NoColor --TerminalIPythonApp.file_to_run='' << eofeofa > log.ipython.nonmag.gw
import test
test.calcGW('pickle', inputpath='inputTestG.py')
launch scattered --ppn 24 --account ss04 --queue batch --walltime=0:30:00
eofeofa
checkRc $? 'ipython nonmag failed'

# Save the execMimic log
mv log log.execMimic.nonmag.gw
checkRc $? 'cannot find nonmag.gw log'
echo testOne: finished nonmag.gw run





















































echo testOne: quit now xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
exit 0

echo ''
echo 'testOne: start magnetic relaxation'
echo ''

# Run magnetic tests
ipython --no-confirm-exit --colors=NoColor --TerminalIPythonApp.file_to_run='' << eofeofb > log.ipython.mag
import test
explore pickle
test.magnetic_wave('pickle', inputpath='inputTestR.py')
launch scattered --ppn 24 --account ss04 --queue short --walltime=0:30:00
eofeofb
checkRc $? 'ipython mag failed'

# Save the execMimic log
mv log log.execMimic.mag
checkRc $? 'cannot find mag log'
echo testOne: finished mag run

# All stderr files must be empty
echo testOne: begin check of stderr
wc -c $(find $workDir -name stderr | sort)
wc -c $(find $workDir -name stderr | sort) | egrep '^0 total$' > /dev/null
checkRc $? 'some stderr file is not empty'
echo testOne: finished check of stderr

# The execMimic logs must have no exc* or err* messages.
echo testOne: start check of execMimic logs
egrep 'exc|error' $workDir/log.execMimic.* > /dev/null
if [[ $? -ne 1 ]]; then
  echo 'log.execMimic.* has exceptions'
  exit 1
fi
echo testOne: finished check of execMimic logs

# Make sure the number of OUTCAR files matches
echo testOne: begin check 1 of OUTCAR files
diff <(cd $refDir; find . -name OUTCAR | sort)  <(cd $workDir; find . -name OUTCAR | sort)
checkRc $? 'OUTCAR files differ'

# Make sure the content of OUTCAR files matches
echo testOne: begin check 2 of OUTCAR files
diff <(cd $refDir; sha1sum $(find . -name OUTCAR | sort))  <(cd $workDir; sha1sum $(find . -name OUTCAR | sort))
checkRc $? 'OUTCAR files differ'
echo testOne: end check of OUTCAR files


# Make sure the number of vasprun.xml files matches
echo testOne: begin check 1 of vasprun.xml files
diff <(cd $refDir; find . -name vasprun.xml | sort)  <(cd $workDir; find . -name vasprun.xml | sort)
checkRc $? 'vasprun.xml files differ'

# Make sure the content of vasprun.xml files matches
echo testOne: begin check 2 of vasprun.xml files
diff <(cd $refDir; sha1sum $(find . -name vasprun.xml | sort))  <(cd $workDir; sha1sum $(find . -name vasprun.xml | sort))
checkRc $? 'vasprun.xml files differ'
echo testOne: end check of vasprun.xml files

echo ''
echo testOne: all done, all OK

