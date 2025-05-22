#!/bin/bash

export lqcd_data_path=/lqcd_data_path

mkdir -p "$lqcd_data_path/Propagators"

#mpirun --allow-run-as-root -n 4 chroma -i ./baryon_2pt.xml
mpirun --allow-run-as-root -n 24 chroma -i ./cpu_inverter.xml
