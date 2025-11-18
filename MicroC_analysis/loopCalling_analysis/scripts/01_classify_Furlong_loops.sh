
mkdir -p commonLoops
mkdir -p undetectableLoops

nLoops=$(ls -1 *PGP*png | wc -l)
echo "Screening ${nLoops} loops"

for n in $(seq 1 1 142);
do
    echo $n

    for resolution in 400 800 1000 2000 4000 8000 10000 20000 40000
    do
	echo $resolution
	for file in $(ls -1 *PGP-${n}_*${resolution}bp.png 2> /dev/null);
	do
	    echo $file
	    open $file
	    
	    read -r -p "What do you wish to do with this file (c - common loop; u - undectable loop)? " answer
	    echo "You typed: $answer"
	    
	    if [[ $answer == "c" ]];
	    then
		mv -v $file commonLoops
		mv -v *PGP-${n}_* commonLoops
		nLoops=$(ls -1 *PGP*png | wc -l)
		echo "${nLoops} loops yet to screen..."
		continue
	    fi
	    if [[ $answer == "u" ]];
	    then
		dir=undetectableLoops
	    fi        
	    echo "Moving ${file} to ${dir}"
	    mv -v ${file} ${dir}
	done
    done
done
