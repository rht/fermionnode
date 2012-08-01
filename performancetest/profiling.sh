set -e
currentdate=$(date "+%H:%M:%S:%m-%d-%Y")
#python -m cProfile -o performancetest/profile-"$currentdate".txt -s cumulative ferminode.py
python -m cProfile -s cumulative ferminode.py > performancetest/profile-"$currentdate".txt 
#python -m cProfile -o performancetest/profile-"$currentdate".txt ferminode.py
