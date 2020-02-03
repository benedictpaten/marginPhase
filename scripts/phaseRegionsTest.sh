
options=${1}

for i in AVG-chr20 AVG-chr7 AVG-chr8 HIGH-chr12 LOW-chr2
do
   time ../scripts/phaseTest.sh ${i} ${options}
done

