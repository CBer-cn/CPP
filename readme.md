# EWSA

## Directory Structure Information

**First category Instances**: Instance sets can be downloaded via https://github.com/hellozhilu/MDMCP.

**Second category Instances**: Instance sets can be downloaded via https://github.com/MMSorensen/CP-Lib.

**EWSA.cpp**: The source code of the EWSA algorithm.

**Supplementary_Material.pdf**: The supplementary material where some illustrations and experiments are presented in forms of figures and tables.

## Running EWSA.cpp

1. Place the EWSA.cpp and data folders in a same directory, and put the instances to be calculated in the data folder.

```
test
|	EWSA.cpp
+----data
|		gauss500-100-1.txt
```

2. Compile EWSA.cpp in Linux environment.

```
g++ -O3 -o EWSA EWSA.cpp
```

3. Run the EWSA .

```
./EWSA data solution
```

**data**: The location of the folder where the test cases are stored, which can be changed as required.

**solution**: The location of the folder where the calculation results are stored, which can be changed as required (the name of results are same as instances).