
#!/bin/bash


minPower=2
maxPower=4

minTrees=50
stepTrees=50
maxTrees=200


mkdir "./R/PI_RESULTS/"

for p in $(seq ${minPower} ${maxPower}); do

for tr in $(seq ${minTrees} ${stepTrees} ${maxTrees}); do

R --vanilla --slave < PI_REG.r --args ${p} ${tr} > result_REG.out
R --vanilla --slave < PI_CLS.r --args ${p} ${tr} > result_CLS.out

done

done









