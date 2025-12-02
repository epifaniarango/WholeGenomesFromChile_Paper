# Scripts for making Figure 1 and related figures

## Filtering Steps and PCA Visualization

First, I merge all the chromosome files into a single file like this and transform it into PLINK format:

```
folder="/your/folder"
output_file="/yourfolder/ModernData.vcf.gz"

vcf_files="${folder}/*.vcf.gz"

# Use --merge to merge by position
bcftools merge -v snps -m2 -M2 --merge all ${vcf_files} -o ${output_file} --force-samples

vcftools --gzvcf ModernData.vcf.gz --plink --out ModernData
```
Here I also keept only biallelic sites.

### Before plotting
We also wanted to crosscheck the call overlap between the SNPChip and the WG sequencing.HO refers to Human Origins (the SNPCHip platform).
Check the R code Comparison.R for more detail.


Then I filter for HWE and the 2 MAF filters (0.01 and 0.05)
```
plink --bfile ModernData --hwe 1e-4  --make-bed --out ModernDataHWE --recode
plink --bfile ModernDataHWE --maf 0.01  --make-bed --out ModernDataHWEmaf01 --recode
plink --bfile ModernDataHWE --maf 0.05  --make-bed --out ModernDataHWEmaf05 --recode
```

Then we filtered to have the same individuals on the SNPChip and WholeGenome(WG) datasets. Before computing PCA, we need to prune to remove linked SNPs. The parameters are different for each dataset:
- --indep-pairwise 50 5 0.4 for WG data
- --indep-pairwise 200 50 0.4 for SNPChip data
  
Pruning is performed like this with any of the parameters before pca:

```
plink --file yourfile --indep-pairwise 200 25 0.4 --out x.tmp
plink --file yourfile --extract x.tmp.prune.in --make-bed --out yourfile.pruned
plink --file yourfile.pruned --pca --out PCA
```

## Fig. S1
For Figure 1, I installed KING with conda and it was pretty straighforward. Use [CodeFigS1.r](https://github.com/epifaniarango/WholeGenomesFromChile_Paper/blob/Figure1andRelated/CodeFigS1.r) to make the plot. 


## ROH 
- I calculated ROH with bcftools using the script [run_roh.sh](https://github.com/epifaniarango/WholeGenomesFromChile_Paper/blob/Figure1andRelated/run_roh.sh) and processed it using the code [Roh.R](https://github.com/epifaniarango/WholeGenomesFromChile_Paper/blob/Figure1andRelated/Roh.R).


## Count variant sites and singletons
With the script [SNPcounts.sh](https://github.com/epifaniarango/WholeGenomesFromChile_Paper/blob/Figure1andRelated/SNPcounts.sh) 

## Code for Figure 1: 
I am quite proud of this [plot](https://github.com/epifaniarango/WholeGenomesFromChile_Paper/blob/Figure1andRelated/CodeFig1.r) 




- To calculate frequencies, I used the flags—-freq and—-het for the F distribution with plink. Check the R script to plot Figure S3. 
- Code for plotting site frequency spectrum on FigureS4 code.
- I calculated Fst with this [code](https://github.com/epifaniarango/Fst_forLargeDatasets), and to plot heatmap and regression, you can check XXX. 




