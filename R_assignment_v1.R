geno<-read_tsv("https://raw.githubusercontent.com/niclizer/BCB_546_UNIX_Assignment/main/assignments/UNIX_Assignment/fang_et_al_genotypes.txt")
view(geno)
snp <- read_tsv("https://raw.githubusercontent.com/niclizer/BCB_546_UNIX_Assignment/main/assignments/UNIX_Assignment/snp_position.txt")
view(snp)
getwd()

# Select SNP_ID, chromosome and Position from SNP_position data set
library(dplyr)
snp_position.selected <- select(snp, c("SNP_ID", "Chromosome", "Position"))
str(snp_position.selected) # 983 rows and 3 columns
glimpse(snp_position.selected) # similar function with str

#next determine the number of rows
nrow(geno)
#2782
ncol(geno)
#986
NROW(na.omit(geno)) 
# 2782
NCOL(na.omit(geno)) 
# 986
geno
object.size(geno)
# 23124352 bytes
object.size(snp)
# 359152 bytes

summary(as.factor(geno$Group))
#TRIPS ZDIPL ZLUXR ZMHUE ZMMIL ZMMLR ZMMMR ZMPBA ZMPIL ZMPJA ZMXCH 
#22    15    17    10   290  1256    27   900    41    34    75 
#ZMXCP ZMXIL ZMXNO ZMXNT ZPERR 
#69     6     7     4     9 
summary(as.factor(snp$Chromosome))
#1=155        10=53        2=127        3=107        4=91        5=122        6=76 #7=97        8=62        9=60 multiple=6  unknown=27 

maize_genotypes <- filter(geno, Group == 'ZMMIL' | Group == 'ZMMLR' | Group == 'ZMMMR')
dim(maize_genotypes)
#1573  986
head(maize_genotypes)
summary(as.factor(maize_genotypes$Group))

teosinte_genotypes <- filter(geno, Group == 'ZMPBA' | Group == 'ZMPIL' | Group == 'ZMPJA')
dim(teosinte_genotypes)
# 975 986
head(teosinte_genotypes)
summary(as.factor(teosinte_genotypes$Group))
#ZMPBA ZMPIL ZMPJA 
#900    41    34


#transpose maize and teosinte
maize_genotypes <- column_to_rownames(maize_genotypes, var = "Sample_ID")
maize_genotypes.tr <- t(maize_genotypes)%>%as.data.frame()%>%rownames_to_column(., var = "SNP_ID")
maize_genotypes.tr <- maize_genotypes.tr[3:nrow(maize_genotypes.tr),]


view(maize_genotypes.tr)
view(snp)

teosinte_genotypes <- column_to_rownames(teosinte_genotypes, var = "Sample_ID")
teosinte_genotypes.tr <- t(teosinte_genotypes)%>%as.data.frame()%>%rownames_to_column(., var = "SNP_ID")
teosinte_genotypes.tr <- teosinte_genotypes.tr[3:nrow(teosinte_genotypes.tr),]


snp_maizegeno <- merge(snp_position.selected, maize_genotypes.tr, by="SNP_ID")
snp_teosintegeno <- merge(snp_position.selected, maize_genotypes.tr, by="SNP_ID")
snp_maizegeno.select <- select(snp_maizegeno, SNP_ID, Chromosome, Position, everything())
snp_teosintegeno.select <- select(snp_teosintegeno, SNP_ID, Chromosome, Position, everything())
view(snp_teosintegeno.select)
dim(snp_maizegeno.select)
#983 1576
summary(as.factor(snp_maizegeno.select$Chromosome))
sum(snp_maizegeno.select$Chromosome == "") 
#0


maize_chrom_inc1 <- subset(snp_maizegeno.select, Chromosome==1)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))
maize_chrom_inc2 <- subset(snp_maizegeno.select, Chromosome==2)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))
maize_chrom_inc3 <- subset(snp_maizegeno.select, Chromosome==3)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))
maize_chrom_inc4 <- subset(snp_maizegeno.select, Chromosome==4)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))
maize_chrom_inc5 <- subset(snp_maizegeno.select, Chromosome==5)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))
maize_chrom_inc6 <- subset(snp_maizegeno.select, Chromosome==6)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))
maize_chrom_inc7 <- subset(snp_maizegeno.select, Chromosome==7)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))
maize_chrom_inc8 <- subset(snp_maizegeno.select, Chromosome==8)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))
maize_chrom_inc9 <- subset(snp_maizegeno.select, Chromosome==9)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))
maize_chrom_inc10 <- subset(snp_maizegeno.select, Chromosome==10)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))

teosinte_chrom_inc1 <- subset(snp_teosintegeno.select, Chromosome==1)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))
teosinte_chrom_inc2 <- subset(snp_teosintegeno.select, Chromosome==2)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))
teosinte_chrom_inc3 <- subset(snp_teosintegeno.select, Chromosome==3)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))
teosinte_chrom_inc4 <- subset(snp_teosintegeno.select, Chromosome==4)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))
teosinte_chrom_inc5 <- subset(snp_teosintegeno.select, Chromosome==5)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))
teosinte_chrom_inc6 <- subset(snp_teosintegeno.select, Chromosome==6)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))
teosinte_chrom_inc7 <- subset(snp_teosintegeno.select, Chromosome==7)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))
teosinte_chrom_inc8 <- subset(snp_teosintegeno.select, Chromosome==8)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))
teosinte_chrom_inc9 <- subset(snp_teosintegeno.select, Chromosome==9)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))
teosinte_chrom_inc10 <- subset(snp_teosintegeno.select, Chromosome==10)%>%arrange(Position)%>%
  mutate_at(4:978, ~replace(., is.na(.), "?"))


maize_chrom_dec1 <- subset(snp_maizegeno.select, Chromosome==1)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
maize_chrom_dec2 <- subset(snp_maizegeno.select, Chromosome==2)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
maize_chrom_dec3 <- subset(snp_maizegeno.select, Chromosome==3)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
maize_chrom_dec4 <- subset(snp_maizegeno.select, Chromosome==4)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
maize_chrom_dec5 <- subset(snp_maizegeno.select, Chromosome==5)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
maize_chrom_dec6 <- subset(snp_maizegeno.select, Chromosome==6)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
maize_chrom_dec7 <- subset(snp_maizegeno.select, Chromosome==7)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
maize_chrom_dec8 <- subset(snp_maizegeno.select, Chromosome==8)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
maize_chrom_dec9 <- subset(snp_maizegeno.select, Chromosome==9)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
maize_chrom_dec10 <- subset(snp_maizegeno.select, Chromosome==10)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))

teosinte_chrom_dec1 <- subset(snp_teosintegeno.select, Chromosome==1)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
teosinte_chrom_dec2 <- subset(snp_teosintegeno.select, Chromosome==2)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
teosinte_chrom_dec3 <- subset(snp_teosintegeno.select, Chromosome==3)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
teosinte_chrom_dec4 <- subset(snp_teosintegeno.select, Chromosome==4)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
teosinte_chrom_dec5 <- subset(snp_teosintegeno.select, Chromosome==5)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
teosinte_chrom_dec6 <- subset(snp_teosintegeno.select, Chromosome==6)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
teosinte_chrom_dec7 <- subset(snp_teosintegeno.select, Chromosome==7)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
teosinte_chrom_dec8 <- subset(snp_teosintegeno.select, Chromosome==8)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
teosinte_chrom_dec9 <- subset(snp_teosintegeno.select, Chromosome==9)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
teosinte_chrom_dec10 <- subset(snp_teosintegeno.select, Chromosome==10)%>%arrange(desc(Position))%>%
  mutate_at(4:978, ~replace(., is.na(.), "-"))
