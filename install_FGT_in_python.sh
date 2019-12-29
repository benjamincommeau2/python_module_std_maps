#!/usr/bin/env python3
#> out.txt
#script -q -a out.txt -c "python3 setup_StdMapPythonC.py build_ext --inplace"
# simply type : bash install_fstoh_in_python.sh : to install python module function chi2
#echo "hello"
#dos2unix out.txt
#script -q -a out.txt -c "echo 'finished'"
#clear
#python3 setup_StdMapPythonC.py build_ext --inplace &>out.txt
#head -n 20 out.txt
if [ $# -eq 0 ]; then
  echo "1" > python_inputs.txt
else
  echo "$1" > python_inputs.txt
fi
python3 setup_StdMapPythonC.py build_ext --inplace 
