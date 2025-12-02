library(data.table)
library(ggplot2)

setwd("~/s3it/0.DataExploration/0.Modern/3.SNP_Comparison/")

map_HO=read.table("GelatoHO_mergedSetAugustBED.map")
map_WH=read.table("../ModernData.map")


overlap <- intersect(map_HO$V2, map_WH$V2)
overlap1=as.data.frame(overlap)

write.table(overlap1, "snps_overlap.txt", row.names = F, col.names = F, quote = F)


famHO=read.table("big_admix_autosomes.fam")
famWG=read.table("WHoverlap.fam")

library(dplyr)

fam=read.table("HOoverlap.fam")
fam2=read.table("WGoverla.fam")

fam2=subset(fam2, fam2$V2 %in%fam$V2)


write.table(famHO, "indOVerHO.txt", row.names = F, col.names = F, quote = F)
write.table(fam2, "indOverWG.txt", row.names = F, col.names = F, quote = F)


mapWG=read.table("WGoverlap1.map")



`%out%` <- function(a,b) ! a %in% b


#as we can see the nomenclature is completely different

mapWG$x=paste(mapWG$V1,mapWG$V4, sep = "_")

mapWG=mapWG[,c(1,5,3,4)]

write.table(mapWG, "WGoverlap1.map", row.names = F, col.names = F, quote = F)
write.table(mapWG$x, "SNPs_nomenclatura.map", row.names = F, col.names = F, quote = F)
mapHO=read.table("HOoverlap.map")
write.table(mapHO$V2, "SNPs_nomenclatura_filterWG.map", row.names = F, col.names = F, quote = F)


write.table(famHO[,1:2], "order.txt", row.names = F, col.names = F, quote = F)

###
mW=read.table("WGmissing.lmiss", header = T)
mH=read.table("HOmissing.lmiss", header = T)

mW=subset(mW, mW$N_MISS == 0)
mH=subset(mH, mH$N_MISS == 0)



###keep the union of snps
mW=subset(mW, mW$SNP %in% mH$SNP)

#create the list to subset with plink 
write.table(mW$SNP, "noMissingSNPlist.txt", row.names = F, col.names = F, quote = F)

#plink --bfile WGoverlap3  --extract noMissingSNPlist.txt --recode --out WGoverlap4 
#plink --bfile HOoverlap1 --extract noMissingSNPlist.txt --recode --out HOoverlap2
###################################
#awk '{$1=$2=$3=$4=$5=$6=""; sub(/^ +/, ""); print}'  WGoverlap3.ped > WGoverlap3.tsv
#awk '{$1=$2=$3=$4=$5=$6=""; sub(/^ +/, ""); print}'  HOoverlap1.ped > HOoverlap1.tsv

#awk '{ for (i = 1; i <= NF; i += 2) { printf "%s%s", $i, $(i + 1); if (i + 2 <= NF) printf " "; } printf "\n"; }' WGoverlap3.tsv > WGoverlap3Merged.tsv
#awk '{ for (i = 1; i <= NF; i += 2) { printf "%s%s", $i, $(i + 1); if (i + 2 <= NF) printf " "; } printf "\n"; }' HOoverlap1.tsv > HOoverlap1Merged.tsv

WGped=read.table("WGoverlap3Merged.tsv")



HOped=read.table("HOoverlap1Merged.tsv")





##############


# Custom function to check if sorted vectors are equal, considering "00" as missing
are_sorted_strings_equal <- function(x, y) {
  x_sorted <- sort(replace(as.character(unlist(strsplit(as.character(x), ''))), x == "00", NA))
  y_sorted <- sort(replace(as.character(unlist(strsplit(as.character(y), ''))), y == "00", NA))
  
  all(is.na(x_sorted) | is.na(y_sorted) | x_sorted == y_sorted)
}

# Create a matrix for the result
result_df <- data.frame(matrix(ncol =ncol(HOped), nrow = nrow(HOped)))
colnames(result_df) <- colnames(HOped)

# Loop through rows and apply the comparison function
for (i in 1:nrow(HOped)) {
  result_df[i, ] <- mapply(are_sorted_strings_equal, HOped[i,], WGped[i,])
}

# View the result

counts_df <- data.frame(
  Row = 1:nrow(result_df),  # Row numbers
  TRUE_Count = rowSums(result_df)  # Count of TRUE values per row
)

write.table(result_df,"MatrixCounts.txt", row.names = F, col.names = F, quote = F)
write.table(counts_df,"CountsJAN.txt", row.names = F, col.names = F, quote = F)



#####
CC=read.table("CountsJan.txt")

fam=read.table("HOoverlap.fam")

counts=cbind(fam[,1:2],CC$V2)
counts$minus=499602-counts$`CC$V2`
write.table(counts,"CountsJANProcessed.txt", row.names = F, col.names = F, quote = F)

# Create a bar plot
ggplot(counts, aes(x = V2, y = minus, fill = V1)) +
  geom_bar(stat = "identity") +
  labs(title = "Bar Plot with V1-based Color", x = "V1", y = "CC.V2 - minus") +
  theme_minimal()
