# GWAS

For 
plink 
### Software:
The tool used for Quality Control and test of association was plink (version 1.90b6.21). While for the study of heritability gcta (version 1.94.1) was used.
All other datamanagement and graph making were done in R

### Quality control

shell prompt
```
plink --bfile GWA-QC --het --allow-no-sex --out GWA-QC 
```

prepartion for reomoval of individuals with a too high heterozygosity rate and NA data.
There are no sorting of missingness, as there are en general a large missingness in the data.
```R
library(tidyverse)
d_het <- read.table("GWA-QC.het",header=T)
d_het$Het<-(d_het$N.NM.-d_het$O.HOM.)/d_het$N.NM.

na_rows<-apply(is.na(d_het),1,any)
rows<-d[na_rows,]
d_het <- na.omit(d_het)
mean_het <- mean(d_het$Het)
sd_het <- sd(d_het$Het)
het_lower <- mean_het - 3*sd_het
het_upper <- mean_het + 3*sd_het
d_m <- d_het %>% filter(Het<het_lower | Het > het_upper)
d_p <- rbind(rows,d_m)
write.table(d_p, file = "wrong_het.txt", col.names = F, row.names = F)
```
removal of individuals with a too high heterozygosity rate or NA data.
shell promps
```bash
plink --bfile GWA-QC --remove wrong_het.txt --allow-no-sex --make-bed --out GWA-QC
```

Identification of duplicated or related indiviudals

shell prompt
```
plink --bfile GWA-QC --indep-pairwise 500kb 5 0.2 --out GWA-QC
```


