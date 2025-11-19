For the samples of the dm6_DPRE2Fab7 condition, you should run the following command to generate the modified dm6 assemblies:

```
for dir in $(ls -1 | grep dm6) ;
do
    echo $dir
    
    bash     ./scripts/01_editFasta.sh &> 01_editFasta.out
    
done
```
