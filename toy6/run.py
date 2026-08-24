#!/usr/bin/env python3
import os
cfile = "toy6_github.cpp"

for tn1 in [0, 1, 4]:

    dirname = f"N6_TAG{tn1}_5"
    os.mkdir(dirname)           # make directory

    cmd = f"cp {cfile} {dirname}/"
    os.system(cmd)              # copy source code

    os.chdir(dirname)
    # compile
    cmd = f"g++ -O3 -std=c++17 -DMMM1={tn1} -o TN={tn1}.out {cfile}"
    ret = os.system(cmd)
    if ret != 0:
        raise RuntimeError(f"Compilation failed for TN1={tn1}")

    # execute
    cmd = f"./TN={tn1}.out > a.dat &"
    os.system(cmd)

    # go back
    os.chdir("..")
