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
