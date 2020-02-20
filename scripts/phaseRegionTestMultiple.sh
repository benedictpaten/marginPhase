
options=${1}

for i in 1 2 3 4 5 6 7 8 9 10
do
   time ../scripts/phaseTest.sh AVG-chr7 ${options} >& out.txt.${i}
done

