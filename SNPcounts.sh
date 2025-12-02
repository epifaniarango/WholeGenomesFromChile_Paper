##########

# 0) Keep autosomal, SNP-only, biallelic sites 
plink --bfile ModernDataHWE \
      --chr 1-22 --snps-only just-acgt --biallelic-only strict \
      --make-bed --out autosome_snps


# export additive genotypes (0/1/2/NA)
plink --bfile autosome_snps --recodeA --out variants_per_indiv

# count SNPs 
awk 'NR==1{next}
{
  c=0; for(i=7;i<=NF;i++) if($i!="NA" && $i>0) c++;
  print $1"\t"$2"\t"c   # FID IID COUNT
}' variants_per_indiv.raw > variants_per_indiv.counts.tsv


#singletons and doubletons, I counted singletons and doubletons but at the end we only use variants and singletons for the final plot


# site frequencies
plink --bfile autosome_snps --freq --out autosome_snps

# make lists; tolerate rounding (0.5≤AC<1.5 is AC≈1; 1.5≤AC<2.5 is AC≈2)
awk 'NR>1{maf=$5; nchr=$6; ac=maf*nchr; if(ac>=0.5 && ac<1.5) print $2}' autosome_snps.frq > singleton.snps
awk 'NR>1{maf=$5; nchr=$6; ac=maf*nchr; if(ac>=1.5 && ac<2.5) print $2}' autosome_snps.frq > doubleton.snps


# singletons
plink --bfile autosome_snps --extract singleton.snps --recodeA --out singletons
awk 'NR==1{next}{c=0; for(i=7;i<=NF;i++) if($i!="NA" && $i>0) c++; print $1"\t"$2"\t"c}' \
  singletons.raw > singletons_per_indiv.tsv

# doubletons
plink --bfile autosome_snps --extract doubleton.snps --recodeA --out doubletons
awk 'NR==1{next}{c=0; for(i=7;i<=NF;i++) if($i!="NA" && $i>0) c++; print $1"\t"$2"\t"c}' \
  doubletons.raw > doubletons_per_indiv.tsv


