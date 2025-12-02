#only mapuche biallelic SNPs 
#I didn't run this script in the end but the commands are here for reference
MAPUCHE_DIR="/Users/ragsdalelab/Documents/PhD/Liftover/Mapuche/data"
OUT="/Users/ragsdalelab/Documents/PhD/Liftover/Mapuche/data"
THREADS=4
MAPUCHE_SAMPLES="/Users/ragsdalelab/Documents/PhD/Liftover/mapuche.samples"

for c in $(seq 1 22); do
  invcf="${MAPUCHE_DIR}/Chile.chr${c}.unphased.vcf.gz"

  bcftools view -S "${MAPUCHE_SAMPLES}" -v snps -m2 -M2 -Oz -o "${OUT}/Mapuche.chr${c}.sites.vcf.gz" "${invcf}" --threads ${THREADS}
  bcftools index -t "${OUT}/Mapuche.chr${c}.sites.vcf.gz"

done




# ------------------------------------------------------
OUT="/Users/ragsdalelab/Documents/PhD/Liftover/Mapuche"
mkdir -p "${OUT}"/{sets,keys,logs}

echo -e "chr\tMapuche_sites" > "${OUT}/counts.tsv"

MAPUCHE_DIR="/Users/ragsdalelab/Documents/PhD/Liftover/Mapuche/data"
for c in $(seq 1 22); do
  invcf="${MAPUCHE_DIR}/Mapuche.chr${c}.sites.vcf.gz"

  bcftools +fill-tags -Ou "${invcf}" -- -t AC,AN,AF \
  | bcftools view -i 'AC>0' -G -Oz -o "${OUT}/sets/Mapuche.chr${c}.sites.vcf.gz" --threads ${THREADS}
  bcftools index -t "${OUT}/sets/Mapuche.chr${c}.sites.vcf.gz"

  
  bcftools query -f '%CHROM:%POS:%REF:%ALT\n' "${OUT}/sets/Mapuche.chr${c}.sites.vcf.gz" \
    | sort -u > "${OUT}/keys/Mapuche.chr${c}.sites"

  nMAP=$(wc -l < "${OUT}/keys/Mapuche.chr${c}.sites" || echo 0)
  printf "chr%s\t%s\n" "${c}" "${nMAP}" >> "${OUT}/counts.tsv"
done

# Totals
totMAP=$(awk 'NR>1 {s+=$2} END{print s+0}' "${OUT}/counts.tsv")
echo -e "TOTAL\t${totMAP}" >> "${OUT}/counts.tsv"



#########COMPARISON WITH 1KGP 

#!/usr/bin/env bash
set -euo pipefail
LC_ALL=C

# --- edit these paths ---
MAP_KEYS="/Users/ragsdalelab/Documents/PhD/Liftover/Mapuche/keys"
KGP_KEYS="/Users/ragsdalelab/Documents/PhD/Liftover/Mapuche37/compare_1kg/perchr"
OUT="/Users/ragsdalelab/Documents/PhD/Liftover/FinalComparison/compare_1kg"
# ------------------------
#I did the keys with a previous script, sorry

mkdir -p "${OUT}"/{perchr,totals}

echo -e "chr\tMapuche\t1KGP_all\tOverlap\tMapuche_only" > "${OUT}/perchr/counts.tsv"



for c in $(seq 1 22); do
  M="${MAP_KEYS}/Mapuche.chr${c}.sites"

  comm -12 "$M" "${KGP_KEYS}/1KGP.ALL.chr${c}.sites" > "${OUT}/perchr/Mapuche_in_1KGP.chr${c}.sites"
  comm -23 "$M" "${KGP_KEYS}/1KGP.ALL.chr${c}.sites" > "${OUT}/perchr/Mapuche_not_in_1KGP.chr${c}.sites"

  nM=$(wc -l < "$M" || echo 0)
  nK=$(wc -l < "${KGP_KEYS}/1KGP.ALL.chr${c}.sites" || echo 0)
  nO=$(wc -l < "${OUT}/perchr/Mapuche_in_1KGP.chr${c}.sites" || echo 0)
  nU=$(wc -l < "${OUT}/perchr/Mapuche_not_in_1KGP.chr${c}.sites" || echo 0)
  printf "chr%s\t%s\t%s\t%s\t%s\n" "$c" "$nM" "$nK" "$nO" "$nU" >> "${OUT}/perchr/counts.tsv"
done

# Totals across chromosomes
cat "${OUT}/perchr/Mapuche_in_1KGP.chr"*.sites       | sort -u > "${OUT}/totals/Mapuche_in_1KGP.all.sites"
cat "${OUT}/perchr/Mapuche_not_in_1KGP.chr"*.sites   | sort -u > "${OUT}/totals/Mapuche_not_in_1KGP.all.sites"
cat "${KGP_KEYS}/1KGP.ALL.chr"*.sites              | sort -u > "${OUT}/totals/1KGP.ALL.all.sites"
cat "${MAP_KEYS}/Mapuche.chr"*.sites                 | sort -u > "${OUT}/totals/Mapuche.all.sites"

Mtot=$(wc -l < "${OUT}/totals/Mapuche.all.sites")
Ktot=$(wc -l < "${OUT}/totals/1KGP.ALL.all.sites")
Otot=$(wc -l < "${OUT}/totals/Mapuche_in_1KGP.all.sites")
Utot=$(wc -l < "${OUT}/totals/Mapuche_not_in_1KGP.all.sites")

printf "\nTOTALS\nMapuche: %'d\n1KGP_all: %'d\nOverlap (Mapuche∩1KGP): %'d (%.2f%% of Mapuche)\nMapuche_only (vs 1KGP): %'d (%.2f%%)\n" \
  "$Mtot" "$Ktot" "$Otot" "$(awk -v a=$Otot -v b=$Mtot 'BEGIN{print 100*a/b}')" \
  "$Utot" "$(awk -v a=$Utot -v b=$Mtot 'BEGIN{print 100*a/b}')"



TOTALS
Mapuche: 6,956,194
1KGP_all: 77,594,210
Overlap (Mapuche∩1KGP): 5,780,896 (83.10% of Mapuche)
Mapuche_only (vs 1KGP): 1,175,298 (16.90%)


################

#!/usr/bin/env bash
set -euo pipefail
LC_ALL=C

OUT="/Users/ragsdalelab/Documents/PhD/FinalComparison/compare_1kg"
KGP_KEYS="/Users/ragsdalelab/Documents/PhD/1KGP/1kg_privNA/keys"

mkdir -p "${OUT}/american_overlap"

echo -e "chr\tMapuche_in_1KGP\tIn_1KGP_NA\tProp_in_NA(%)" > "${OUT}/american_overlap/counts.tsv"

for c in $(seq 1 22); do
  map_in_all="${OUT}/perchr/Mapuche_in_1KGP.chr${c}.sites"
  kgp_na="${KGP_KEYS}/NA.chr${c}.sites"

  [[ -f "$map_in_all" && -f "$kgp_na" ]] || { echo "Missing chr${c} files"; continue; }

  comm -12 "$map_in_all" "$kgp_na" > "${OUT}/american_overlap/Mapuche_in_1KGP_NA.chr${c}.sites"

  n_all=$(wc -l < "$map_in_all" || echo 0)
  n_na=$(wc -l < "${OUT}/american_overlap/Mapuche_in_1KGP_NA.chr${c}.sites" || echo 0)

  prop=$(awk -v a=$n_na -v b=$n_all 'BEGIN{if(b>0) printf("%.2f",100*a/b); else print "0"}')
  printf "chr%s\t%s\t%s\t%s\n" "$c" "$n_all" "$n_na" "$prop" >> "${OUT}/american_overlap/counts.tsv"
done

# Combine across chromosomes
cat "${OUT}/american_overlap/Mapuche_in_1KGP_NA.chr"*.sites | sort -u > "${OUT}/american_overlap/Mapuche_in_1KGP_NA.all.sites"

n_all_total=$(cat "${OUT}/perchr/Mapuche_in_1KGP.chr"*.sites | sort -u | wc -l)
n_na_total=$(wc -l < "${OUT}/american_overlap/Mapuche_in_1KGP_NA.all.sites")
prop_total=$(awk -v a=$n_na_total -v b=$n_all_total 'BEGIN{printf("%.2f",100*a/b)}')

echo
echo "TOTAL: Mapuche∩1KGP = $n_all_total sites"
echo "       Of these, in 1KGP American (NA) = $n_na_total sites ($prop_total%)"



TOTAL: Mapuche∩1KGP =  5780896 sites
       Of these, in 1KGP American (NA) =  5713255 sites (98.83%)




##########Of the non-shared sites, how many are singletons in Mapuche?


UNI="/Users/ragsdalelab/Documents/PhD/Liftover/MapucheExcluded/compare_1kg/totals/Mapuche_not_in_1KGP.all.sites"
#!/usr/bin/env bash
set -euo pipefail
LC_ALL=C

# --- EDIT THESE PATHS ---
MAP_VCF_DIR="/Users/ragsdalelab/Documents/PhD/Liftover/Mapuche/data"   
UNIQUE_KEYS="/Users/ragsdalelab/Documents/PhD/Liftover/FinalComparison/compare_1kg/totals/Mapuche_not_in_1KGP.all.sites"  
OUT="/Users/ragsdalelab/Documents/PhD/Liftover/FinalComparison/compare_1kg/unique_vs1kg_vcf"
THREADS=4
# ------------------------

mkdir -p "${OUT}/perchr" "${OUT}/logs" "${UNI_DIR}/totals"

conda activate bcftools

# 1) Build 1-based CHROM\tPOS targets from keys (positions only)
POS="${OUT}/unique_vs1kg.targets.pos"
awk -F':' 'BEGIN{OFS="\t"} {print $1,$2}' "${UNIQUE_KEYS}" \
 | sort -k1,1V -k2,2n > "${POS}"

# 2) Extract those sites from each chromosome VCF, fill AC/AN/AF, write per-chr VCFs
for c in $(seq 1 22); do
  in_vcf="${MAP_VCF_DIR}/Mapuche.chr${c}.sites.vcf.gz"
  out_vcf="${OUT}/perchr/Mapuche.unique_vs1kg.chr${c}.vcf.gz"
  [[ -f "$in_vcf" ]] || { echo "Missing $in_vcf" >&2; continue; }

  echo "[chr${c}] subsetting to Mapuche-unique (vs 1KGP)…"
  bcftools view -T "${POS}" -Ou "$in_vcf" \
  | bcftools +fill-tags -Ou -- -t AC,AN,AF \
  | bcftools view -Oz -o "$out_vcf" --threads ${THREADS}
  bcftools index -f "$out_vcf"
done

# 3) Concatenate all chromosomes → one VCF of Mapuche-unique (vs 1KGP)
ALL_VCF="${OUT}/Mapuche.unique_vs1kg.ALL.vcf.gz"
bcftools concat -Oz --threads ${THREADS} \
  -o "${ALL_VCF}" "${OUT}/perchr"/Mapuche.unique_vs1kg.chr*.vcf.gz
bcftools index -f "${ALL_VCF}"

# 4) Allele-count profile (singletons/doubletons/MAC bins)
total=$(bcftools view -H "${ALL_VCF}" | wc -l | awk '{print $1}')
single=$(bcftools view -i 'INFO/AC==1' -H "${ALL_VCF}" | wc -l | awk '{print $1}')
double=$(bcftools view -i 'INFO/AC==2' -H "${ALL_VCF}" | wc -l | awk '{print $1}')
mac_le2=$(bcftools view -i 'INFO/AC<=2' -H "${ALL_VCF}" | wc -l | awk '{print $1}')
mac_le5=$(bcftools view -i 'INFO/AC<=5' -H "${ALL_VCF}" | wc -l | awk '{print $1}')

pct() { awk -v a="$1" -v b="$2" 'BEGIN{if(b>0) printf("%.2f",100*a/b); else print "0.00"}'; }

echo
echo "Mapuche-unique (vs 1KGP) VCF: ${ALL_VCF}"
printf "Total unique sites:   %'d\n" "$total"
printf "Singletons (AC=1):    %'d (%s%%)\n" "$single"  "$(pct "$single" "$total")"
printf "Doubletons (AC=2):    %'d (%s%%)\n" "$double"  "$(pct "$double" "$total")"
printf "MAC ≤ 2:              %'d (%s%%)\n" "$mac_le2" "$(pct "$mac_le2" "$total")"
printf "MAC ≤ 5:              %'d (%s%%)\n" "$mac_le5" "$(pct "$mac_le5" "$total")"


Total unique sites:   1,175,298
Singletons (AC=1):    907,400 (77.21%)
Doubletons (AC=2):    145,337 (12.37%)
MAC ≤ 2:              1,052,737 (89.57%)
MAC ≤ 5:              1,129,562 (96.11%)
