# GWAS
## Software:
The tool used for Quality Control and test of association was plink (version 1.90b6.21). While for the study of heritability gcta (version 1.94.1) was used.
All other datamanagement and graph making were done in R

## Quality control

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

```bash
plink --bfile GWA-QC --remove wrong_het.txt --allow-no-sex  --make-bed --out GWA-QC-het
```

### Identification of duplicated or related indiviudals

```
 plink --allow-no-sex --bfile GWA-QC-het --indep-pairwise 500kb 5 0.2 --out GWA-QC-het

plink --bfile GWA-QC-het --extract GWA-QC-het.prune.in --genome --min 0.185 --allow-no-sex --out GWA-QC-het
```
Identification of the IBD individuals in R
```R
ibd <- read.table('GWA-QC-het.genome', header = TRUE)
members <- ibd$FID1
members <- unique(members)
write.table(cbind(members,members), file = 'wrong_ibd.txt', col.names = F, row.names = F)
```
The amount of individuals that were IBD were foun with
```
length(members)
```
Removal of the IBD individuals
```bash
plink --bfile  GWA-QC-het --remove wrong_ibd.txt --allow-no-sex  --make-bed --out GWA-QC-ibd
```
## Addition of phenotype
Prepartion of file with individuals with no selfreported height
```R
library(tidyverse)
height<-read.table('height.txt')
indi <- read.table("GWA-QC-ibd.fam",header=F)
height_out<-indi%>% anti_join(height,by=c("V1"="V1"))
write.table(height_out, file = "height_out.txt", col.names = F, row.names = F)
```
removal of individuals with no selfreported height
```bash
 plink --bfile  GWA-QC-ibd --remove height_out.txt --allow-no-sex  --make-bed --out GWA-height
```
Addition of phenotype to a new fam file
```R
fam_data<-read.table("GWA-height.fam",header=F)
fam_data<-merge(fam_data,height,by="V1")
fam_data$V6<-fam_data$V2.y
fam_data<-fam_data[,!(names(fam_data)%in% "V2.y")]
colnames(fam_data) <- c("FID", "IID", "PID", "MID", "Sex", "Phenotype")
write.table(fam_data,"GWA-height_cor.fam" , sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
```
Adding the phenotype to the apropriate dataformat
```bash
plink --bfile GWA-height --pheno GWA-height_cor.fam --pheno-name Phenotype --allow-no-sex  --make-bed --out GWA-full-height
```

## PCA
Creation of the principle components
shell prompt
```
plink --allow-no-sex --bfile GWA-QC-ibd --indep-pairwise 500kb 5 0.2 --out GWA-QC-ibd

plink --allow-no-sex --bfile GWA-QC-ibd --extract GWA-QC-ibd.in --pca 20 --out GWA-QC-ibd
```
```R
library(tidyverse)
eigenvec_height<-read.table('GWA-full-height.eigenvec',head=F)
info<-read.table('GWA-full-height.fam', head=F)
eigenvec_height<-merge(eigenvec_height,info,by="V1")
sex_mapping <- c("Ambigous", "Male", "Female")
eigenvec_height$Sex<-sex_mapping[eigenvec_height$V5.y+1]
eigenvec_height$Height<-eigenvec_height$V6.y
eigenvec_height$PC1<-eigenvec_height$V3.x
eigenvec_height$PC2<-eigenvec_height$V4.x
ggplot(eigenvec_height,aes(x=PC1,y=PC2,color=Sex))+
  geom_point()+
  scale_color_manual(values = c("black","blue", "red"))
ggplot(eigenvec_height,aes(x=PC1,y=PC2,color=Height))+
  geom_point()
```
## Association study
```bash
plink --bfile GWA-full-height --allow-no-sex --linear  --out GWA-full-height
```

```R
library(tidyverse)
library("qqman")
vignette("qqman")
assoc <- read.table('GWA-full-height.assoc.linear', head=T)
assoc <- na.omit(assoc) %>% subset(CHR<=23)
manhattan(assoc,suggestiveline=F)
qq(assoc$P)
```

The amount of associated SNPs were found with the following
```
snp<- assoc%>% subset(P<0.00000005)
nrow(snp)
snp_bonf<- assoc%>% subset(P<(0.00000005/nrow(assoc)))
nrow(snp_bonf)
```
## Heritability
```bash
plink --make-grm-gz --bfile GWA-full-height --allow-no-sex  --out full-height
gcta64 --grm-gz full-height --pheno heights.phen --reml --out heights_test
```

## Ploygenic scores prediction
Preparation of discovery and target samples 
```R
library(data.table)
library(tidyverse
fam <- fread("GWA-full-height.fam", header = FALSE)
colnames(fam) <- c("FID", "IID", "PID", "MID", "Sex", "Phenotype")
set.seed(123) 
shuffled_indices <- sample(nrow(fam))
num_discovery <- round(0.6 * nrow(fam))
discovery_indices <- shuffled_indices[1:num_discovery]
target_indices <- shuffled_indices[(num_discovery + 1):nrow(fam)]
discovery_fam <- fam[discovery_indices, .(FID, IID)]
target_fam <- fam[target_indices, .(FID, IID)]
fwrite(discovery_fam, "discovery_ids.txt", quote = FALSE, sep = "\t", col.names = FALSE)
fwrite(target_fam, "target_ids.txt", quote = FALSE, sep = "\t", col.names = FALSE)
```
Creation of bedfiles for the discovery and target samples
```bash
plink --bfile GWA-full-height --keep discovery_ids.txt --allow-no-sex  --make-bed --out discovery_sample
plink --bfile GWA-full-height --keep target_ids.txt --allow-no-sex --make-bed --out target_sample
```
Linear association of the discovery sample to get weights for the SNPs
```
plink --bfile  discovery_sample -assoc --allow-no-sex --linear  --out gwas_discovery_height
```
The weight of the specific SNPs are sorted into their own file for further calculation of the polygenic score. For SNPs with more than one possible base have their weight averaged for the indiviudal SNP
```R
assoc <- read.table('gwas_discovery_height.assoc.linear', head=T)
assoc <- na.omit(assoc) %>% subset(CHR<=23)
snp<- assoc%>% subset(P<0.00000005)
score_data<-snp%>%select(SNP,A1,BETA)
score_data <- score_data %>%
  group_by(SNP, A1) %>%
  summarise(BETA = mean(BETA)) %>%
  ungroup()
write.table(score_data,"snps_weights.txt",row.names=F,quote=F,sep="\t")
```
Calculation of polygenic score of the target sample
```bash
plink --bfile target_sample --score snp_weights.txt 1 2 3 sum --out pgs_results_height
```
calculation of the prediction accuracy of the polygenic score.
```R
pgs <- fread("pgs_results_height.profile")
phenotypes <- fread("GWA-height_cor.fam")
data <- merge(pgs, phenotypes, by = "FID")
model <- lm(Phenotype ~ SCORESUM, data = data)
summary(model)
```

