#! /bin/bash
# nohup bash allbatch.sh > ./allbatch.txt 2>&1 &

cd single_pn
echo single_pn
Rscript data.R --s=25 > ./log_25.txt 2>&1
# Rscript data.R --s=50 > ./log_50.txt 2>&1
# Rscript data.R --s=100 > ./log_100.txt 2>&1
#Rscript data.R --s=200 > ./log_200.txt 2>&1


cd ../double_pn
echo double_pn
Rscript data.R --s=25 > ./log_25.txt 2>&1
# Rscript data.R --s=50 > ./log_50.txt 2>&1
# Rscript data.R --s=100 > ./log_100.txt 2>&1
#Rscript data.R --s=200 > ./log_200.txt 2>&1


cd ../single_bn
echo single_bn
Rscript data.R --s=25 > ./log_25.txt 2>&1
# Rscript data.R --s=50 > ./log_50.txt 2>&1
# Rscript data.R --s=100 > ./log_100.txt 2>&1
#Rscript data.R --s=200 > ./log_200.txt 2>&1

cd ../double_bn
echo double_bn
Rscript data.R --s=25 > ./log_25.txt 2>&1
# Rscript data.R --s=50 > ./log_50.txt 2>&1
# Rscript data.R --s=100 > ./log_100.txt 2>&1
#Rscript data.R --s=200 > ./log_200.txt 2>&1

cd ../double_hn
echo double_hn
Rscript data.R --s=25 > ./log_25.txt 2>&1
# Rscript data.R --s=50 > ./log_50.txt 2>&1
# Rscript data.R --s=100 > ./log_100.txt 2>&1
#Rscript data.R --s=200 > ./log_200.txt 2>&1

cd ../

echo END
