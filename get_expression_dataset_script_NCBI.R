# we have dataframes (SraRunTable) with Run information and Sample.Name. 
# for analysis select a subset of 50 runs for each dataset
GSE156451 <- sample(GSE156451_SraRunTable$Run,50)
GSE185055 <- sample(GSE185055_SraRunTable$Run,50)

# save the Run information as csv files
write.csv(GSE156451,"GSE156451.csv",row.names = FALSE,col.names = FALSE)
write.csv(GSE185055,"GSE185055.csv",row.names = FALSE,col.names = FALSE)

# further steps use the SRA toolkit from NCBI using the command line.
# prefetch
#!/usr/bin/bash
echo Hello
for entry in `cat GSE156451.csv`
do
prefetch -p $entry
done
echo done 
# fasterq-dump
#!/usr/bin/bash
counter=0
for file in $(ls)
do
if [ -d $file ]
then
echo $file is directory
fasterq-dump --concatenate-reads $file -e 2 -b 100MB -c 1000MB -m 2000MB -p
elif [ -f $file ]
then     
echo $file is file        
((counter--));
fi
sleep 5
((counter++));
done
free -m | awk 'NR==2{printf "Memory Usage: %s/%sMB (%.2f)", $3,$2,$3*100/$2}'
echo \n
echo $counter directories mapped
