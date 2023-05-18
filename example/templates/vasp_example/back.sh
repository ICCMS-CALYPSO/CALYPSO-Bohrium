#!/bin/bash

for i in {1..9}
do
	cp data/step002.pop00${i}/POSCAR POSCAR_${i}
	cp data/step002.pop00${i}/CONTCAR CONTCAR_${i}
	cp data/step002.pop00${i}/OUTCAR OUTCAR_${i}
done

for i in {10..50}
do
	cp data/step003.pop0${i}/POSCAR POSCAR_${i}
	cp data/step003.pop0${i}/CONTCAR CONTCAR_${i}
	cp data/step003.pop0${i}/OUTCAR OUTCAR_${i}
done
