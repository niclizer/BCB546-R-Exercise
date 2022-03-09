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

