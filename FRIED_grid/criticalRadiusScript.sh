#!/bin/bash



cd ..
lastfile="$(ls -Art radial*.dat | tail -n 1)"
lastfile="radial0018.dat"

oldno="2.22507-308"
newno="1.e-30"

echo "last file was " ${lastfile} ${lastfile2}

sed -i -e "s/${oldno}/${newno}/g" ${lastfile}
echo "done first sed"

dummy="inputfile"

outfile="plotCritRadiusRun.py"

cd critPlotPackage
cp plotCritRadius.py plotCritRadiusRun.py

echo "doing second sed", ${dummy}, ${lastfile}
sed -i -e "s/${dummy}/${lastfile}/g" ${outfile}


#"plotCritRadiusRun.py"

python plotCritRadiusRun.py
