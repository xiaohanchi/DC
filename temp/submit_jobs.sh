#!/bin/bash

n_jobs=148
# my_array=(15 16 18 19 223 224 226 227)

# for i in "${my_array[@]}"; do
for i in $(seq 1 ${n_jobs}); do

	mkdir -p jobtype$i/results

	cp main.R main_tmp.R
	perl -pi -e "s/sc00/${i}/g" main_tmp.R
	mv main_tmp.R ./jobtype$i/main.R


	cp runjobs.lsf runjobs_tmp.lsf
	perl -pi -e "s/folder00/jobtype$i/g; s/idx00/$i/g" runjobs_tmp.lsf
	if test $i -gt $((n_jobs / 3)); then
		perl -pi -e "s/short/medium/g; s/3:00/4:00/g" runjobs_tmp.lsf
	fi
	mv runjobs_tmp.lsf ./jobtype$i/runjobs.lsf

done

echo "DONE Creating Files."

# for i in "${my_array[@]}"; do
for i in $(seq 1 ${n_jobs}); do

	bsub <./jobtype$i/runjobs.lsf 

done

echo "DONE Submitting Jobs."
