
#Title
library(dplyr)
library(tidyverse)
library(ggplot2)
library(knitr)

geno<-read_tsv("https://raw.githubusercontent.com/niclizer/BCB_546_UNIX_Assignment/main/assignments/UNIX_Assignment/fang_et_al_genotypes.txt")
snp <- read_tsv("https://raw.githubusercontent.com/niclizer/BCB_546_UNIX_Assignment/main/assignments/UNIX_Assignment/snp_position.txt")

# Select SNP_ID, chromosome and Position from SNP_position data set
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
#ZMMIL ZMMLR ZMMMR 
#290  1256    27

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
view(geno)

teosinte_genotypes <- column_to_rownames(teosinte_genotypes, var = "Sample_ID")
teosinte_genotypes.tr <- t(teosinte_genotypes)%>%as.data.frame()%>%rownames_to_column(., var = "SNP_ID")
teosinte_genotypes.tr <- teosinte_genotypes.tr[3:nrow(teosinte_genotypes.tr),]

view(snp_maizegeno.select)
snp_maizegeno <- merge(snp_position.selected, maize_genotypes.tr, by="SNP_ID")
snp_teosintegeno <- merge(snp_position.selected, teosinte_genotypes.tr, by="SNP_ID")
snp_maizegeno.select <- select(snp_maizegeno, SNP_ID, Chromosome, Position, everything())
snp_teosintegeno.select <- select(snp_teosintegeno, SNP_ID, Chromosome, Position, everything())



maize_chrom_inc1 <- subset(snp_maizegeno.select, Chromosome==1)%>%arrange(as.numeric(Position))
maize_chrom_inc2 <- subset(snp_maizegeno.select, Chromosome==2)%>%arrange(as.numeric(Position))
maize_chrom_inc3 <- subset(snp_maizegeno.select, Chromosome==3)%>%arrange(as.numeric(Position))
maize_chrom_inc4 <- subset(snp_maizegeno.select, Chromosome==4)%>%arrange(as.numeric(Position))
maize_chrom_inc5 <- subset(snp_maizegeno.select, Chromosome==5)%>%arrange(as.numeric(Position))
maize_chrom_inc6 <- subset(snp_maizegeno.select, Chromosome==6)%>%arrange(as.numeric(Position))
maize_chrom_inc7 <- subset(snp_maizegeno.select, Chromosome==7)%>%arrange(as.numeric(Position))
maize_chrom_inc8 <- subset(snp_maizegeno.select, Chromosome==8)%>%arrange(as.numeric(Position))
maize_chrom_inc9 <- subset(snp_maizegeno.select, Chromosome==9)%>%arrange(as.numeric(Position))
maize_chrom_inc10 <- subset(snp_maizegeno.select, Chromosome==10)%>%arrange(as.numeric(Position))

teosinte_chrom_inc1 <- subset(snp_teosintegeno.select, Chromosome==1)%>%arrange(as.numeric(Position))
teosinte_chrom_inc2 <- subset(snp_teosintegeno.select, Chromosome==2)%>%arrange(as.numeric(Position))
teosinte_chrom_inc3 <- subset(snp_teosintegeno.select, Chromosome==3)%>%arrange(as.numeric(Position))
teosinte_chrom_inc4 <- subset(snp_teosintegeno.select, Chromosome==4)%>%arrange(as.numeric(Position))
teosinte_chrom_inc5 <- subset(snp_teosintegeno.select, Chromosome==5)%>%arrange(as.numeric(Position))
teosinte_chrom_inc6 <- subset(snp_teosintegeno.select, Chromosome==6)%>%arrange(as.numeric(Position))
teosinte_chrom_inc7 <- subset(snp_teosintegeno.select, Chromosome==7)%>%arrange(as.numeric(Position))
teosinte_chrom_inc8 <- subset(snp_teosintegeno.select, Chromosome==8)%>%arrange(as.numeric(Position))
teosinte_chrom_inc9 <- subset(snp_teosintegeno.select, Chromosome==9)%>%arrange(as.numeric(Position))
teosinte_chrom_inc10 <- subset(snp_teosintegeno.select, Chromosome==10)%>%arrange(as.numeric(Position))


write.csv(maize_chrom_inc1,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//maize_chrom_inc1")
write.csv(maize_chrom_inc2,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//maize_chrom_inc2")
write.csv(maize_chrom_inc3,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//maize_chrom_inc3")
write.csv(maize_chrom_inc4,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//maize_chrom_inc4")
write.csv(maize_chrom_inc5,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//maize_chrom_inc5")
write.csv(maize_chrom_inc6,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//maize_chrom_inc6")
write.csv(maize_chrom_inc7,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//maize_chrom_inc7")
write.csv(maize_chrom_inc8,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//maize_chrom_inc8")
write.csv(maize_chrom_inc9,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//maize_chrom_inc9")
write.csv(maize_chrom_inc10,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//maize_chrom_inc10")
write.csv(teosinte_chrom_inc1,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//teosinte_chrom_inc1")
write.csv(teosinte_chrom_inc2,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//teosinte_chrom_inc2")
write.csv(teosinte_chrom_inc3,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//teosinte_chrom_inc3")
write.csv(teosinte_chrom_inc4,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//teosinte_chrom_inc4")
write.csv(teosinte_chrom_inc5,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//teosinte_chrom_inc5")
write.csv(teosinte_chrom_inc6,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//teosinte_chrom_inc6")
write.csv(teosinte_chrom_inc7,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//teosinte_chrom_inc7")
write.csv(teosinte_chrom_inc8,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//teosinte_chrom_inc8")
write.csv(teosinte_chrom_inc9,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//teosinte_chrom_inc9")
write.csv(teosinte_chrom_inc10,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//teosinte_chrom_inc10")







data.frame(lapply(data, gsub, pattern = "[?]", replacement = "-"))

view(maize_chom_dec1_1)

maize_replaced<- data.frame(lapply(snp_maizegeno.select, gsub, pattern = "[?]", replacement = "-"))
teosente_replaced <- data.frame(lapply(snp_teosintegeno.select, gsub, pattern = "[?]", replacement = "-"))

maize_chrom_dec1 <- subset(maize_replaced, Chromosome==1)%>%arrange(desc(as.numeric(Position)))
maize_chrom_dec2 <- subset(maize_replaced, Chromosome==2)%>%arrange(desc(as.numeric(Position)))
maize_chrom_dec3 <- subset(maize_replaced, Chromosome==3)%>%arrange(desc(as.numeric(Position)))
maize_chrom_dec4 <- subset(maize_replaced, Chromosome==4)%>%arrange(desc(as.numeric(Position)))
maize_chrom_dec5 <- subset(maize_replaced, Chromosome==5)%>%arrange(desc(as.numeric(Position)))
maize_chrom_dec6 <- subset(maize_replaced, Chromosome==6)%>%arrange(desc(as.numeric(Position)))
maize_chrom_dec7 <- subset(maize_replaced, Chromosome==7)%>%arrange(desc(as.numeric(Position)))
maize_chrom_dec8 <- subset(maize_replaced, Chromosome==8)%>%arrange(desc(as.numeric(Position)))
maize_chrom_dec9 <- subset(maize_replaced, Chromosome==9)%>%arrange(desc(as.numeric(Position)))
maize_chrom_dec10 <- subset(maize_replaced, Chromosome==10)%>%arrange(desc(as.numeric(Position)))

teosinte_chrom_dec1 <- subset(teosente_replaced, Chromosome==1)%>%arrange(desc(as.numeric(Position)))
teosinte_chrom_dec2 <- subset(teosente_replaced, Chromosome==2)%>%arrange(desc(as.numeric(Position)))
teosinte_chrom_dec3 <- subset(teosente_replaced, Chromosome==3)%>%arrange(desc(as.numeric(Position)))
teosinte_chrom_dec4 <- subset(teosente_replaced, Chromosome==4)%>%arrange(desc(as.numeric(Position)))
teosinte_chrom_dec5 <- subset(teosente_replaced, Chromosome==5)%>%arrange(desc(as.numeric(Position)))
teosinte_chrom_dec6 <- subset(teosente_replaced, Chromosome==6)%>%arrange(desc(as.numeric(Position)))
teosinte_chrom_dec7 <- subset(teosente_replaced, Chromosome==7)%>%arrange(desc(as.numeric(Position)))
teosinte_chrom_dec8 <- subset(teosente_replaced, Chromosome==8)%>%arrange(desc(as.numeric(Position)))
teosinte_chrom_dec9 <- subset(teosente_replaced, Chromosome==9)%>%arrange(desc(as.numeric(Position)))
teosinte_chrom_dec10 <- subset(teosente_replaced, Chromosome==10)%>%arrange(desc(as.numeric(Position)))

write.csv(maize_chrom_dec1,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//maize_chrom_dec1")
write.csv(maize_chrom_dec2,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//maize_chrom_dec2")
write.csv(maize_chrom_dec3,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//maize_chrom_dec3")
write.csv(maize_chrom_dec4,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//maize_chrom_dec4")
write.csv(maize_chrom_dec5,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//maize_chrom_dec5")
write.csv(maize_chrom_dec6,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//maize_chrom_dec6")
write.csv(maize_chrom_dec7,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//maize_chrom_dec7")
write.csv(maize_chrom_dec8,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//maize_chrom_dec8")
write.csv(maize_chrom_dec9,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//maize_chrom_dec9")
write.csv(maize_chrom_dec10,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//maize_chrom_dec10")

write.csv(teosinte_chrom_dec1,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//teosinte_chrom_dec1")
write.csv(teosinte_chrom_dec2,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//teosinte_chrom_dec2")
write.csv(teosinte_chrom_dec3,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//teosinte_chrom_dec3")
write.csv(teosinte_chrom_dec4,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//teosinte_chrom_dec4")
write.csv(teosinte_chrom_dec5,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//teosinte_chrom_dec5")
write.csv(teosinte_chrom_dec6,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//teosinte_chrom_dec6")
write.csv(teosinte_chrom_dec7,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//teosinte_chrom_dec7")
write.csv(teosinte_chrom_dec8,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//teosinte_chrom_dec8")
write.csv(teosinte_chrom_dec9,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//teosinte_chrom_dec9")
write.csv(teosinte_chrom_dec10,"C:/Users/nel0208144/Documents/GitHub/EEOB546_R_lesson//teosinte_chrom_dec10")



#maize_snp_bar <- ggplot(data = snp_maizegeno.select) + geom_bar(mapping=aes(x=Chromosome)) + 
 # ggtitle("SNP distribution in Chromosome - Maize genotypes")
#maize_snp_bar + theme(plot.title = element_text(color = "blue", size = 12, face = "bold", hjust = 0.5))
#teosinte_snp_bar <- ggplot(data = snp_teosintegeno.select) + geom_point(mapping=aes(x=Chromosome, y=Position)) +
 # ggtitle("SNP distribution in Chromosome - Teosinte genotypes")
#teosinte_snp_bar + theme(plot.title = element_text(color = "blue", size = 12, face = "bold", hjust = 0.5))

maize_snp_bar <- ggplot(data = snp_maizegeno.select) + geom_bar(mapping=aes(x=Chromosome)) + 
  ggtitle("Miaze SNP distribution across Chromosomes")
maize_snp_bar + theme(plot.title = element_text(color = "blue", size = 12, face = "bold", hjust = 0.5))

view(geno)
view(snp_teosintegeno.select)

teosinte_snp_bar <- ggplot(data = snp_teosintegeno.select) + geom_bar(mapping=aes(x=Chromosome)) + 
  ggtitle("Teosinte SNP distribution across Chromosomes")
teosinte_snp_bar + theme(plot.title = element_text(color = "blue", size = 12, face = "bold", hjust = 0.5))

maize_snp_point <- ggplot(data = snp_maizegeno.select) + geom_point(mapping=aes(x=Chromosome, y=Position)) + 
  ggtitle("Miaze SNP distribution across Chromosomes")
maize_snp_point + theme(plot.title = element_text(color = "blue", size = 12, face = "bold", hjust = 0.5))

teosinte_snp_point <- ggplot(data = snp_teosintegeno.select) + geom_point(mapping=aes(x=Chromosome, y=Position)) + 
  ggtitle("Teosinte SNP distribution across Chromosomes")
teosinte_snp_point + theme(plot.title = element_text(color = "blue", size = 12, face = "bold", hjust = 0.5))

summary(as.factor(snp_maizegeno.select$Chromosome))
summary(as.factor(snp_teosintegeno.select$Chromosome))

#sapply(as.numeric(snp_teosintegeno.select$Position), range)
#sapply(snp_teosintegeno.select$Position, max)


###########
snp2 <- snp[c(1,3,4,2,5:15)]
snp_position.selected
snp_maizegeno
snp_teosintegeno
testplot <- ggplot(data = snp_position.selected[(as.numeric(snp$Chromosome))]) +
  geom_bar(mapping = aes(as.numeric(Chromosome),fill=Chromosome)) + 
  scale_x_discrete(limit=c(1:10)) + labs(x="chromosome number", y="position")
print(testplot)
ggplot(data = snp_position.selected[!is.na(as.numeric(snp$Chromosome)),]) +   geom_bar(mapping = aes(x = as.numeric(Chromosome), fill=Chromosome)) + scale_x_discrete(limit=c(1:10))+ labs(x = "Chromosome number", y="No. of polymorphism position")
ggplot(data = snp2[!is.na(as.numeric(snp$Chromosome)),]) +   geom_bar(mapping = aes(x = as.numeric(Chromosome), fill=Chromosome))



#dead section below
maize_genotypes <- filter(geno, Group == 'ZMMIL' | Group == 'ZMMLR' | Group == 'ZMMMR')
teosinte_genotypes <- filter(geno, Group == 'ZMPBA' | Group == 'ZMPIL' | Group == 'ZMPJA')
maize_pivot <- maize_genotypes %>% pivot_longer(!c(Sample_ID, JG_OTU, Group), names_to="SNP_ID", values_to= "NT")
teosinte_pivot <- teosinte_genotypes %>% pivot_longer(!c(Sample_ID, JG_OTU, Group), names_to="SNP_ID", values_to= "NT")
snp_maize_pivot <- merge(maize_pivot, snp_position.selected, by="SNP_ID")
snp_teosinte_pivot <- merge(teosinte_pivot, snp_position.selected, by="SNP_ID")
num_snp_maize_pivot <- (snp_maize_pivot[!is.na(as.numeric(snp_maize_pivot$Chromosome)),])
num_snp_teosinte_pivot <- (snp_teosinte_pivot[!is.na(as.numeric(snp_teosinte_pivot$Chromosome)),])


maize_identify <- snp_maizegeno.select%>% select(SNP_ID, Chromosome, Position) %>% mutate(Species = "Maize")
teosinte_identify <- snp_teosintegeno.select%>% select(SNP_ID, Chromosome, Position) %>% mutate(Species = "Teosinte")
genobound <- bind_rows(maize_identify, teosinte_identify)
view(genobound)

fang_pivot_longer <- geno %>% pivot_longer(!c(Sample_ID, JG_OTU, Group), names_to="SNP_ID", values_to= "NT")
snp_fang <- merge(fang_pivot_longer, snp, by="SNP_ID")
num_snp_fang <- (snp_fang[!is.na(as.numeric(snp_fang$Chromosome)),])

SNPs_per_chrom <- (genobound %>% select(SNP_ID, Chromosome, Position) %>% drop_na() %>% ggplot()+
      geom_bar(mapping = aes(as.numeric(Chromosome)), fill = Species, color = Species) +
      labs(x = "Chromosome", y = "Total SNPs") + ggtitle("SNPs in Chromosome Position") +
      theme(plot.title = element_text(hjust = 0.5)))


SNPs_per_chrom3 <- ggplot(data = genobound) + geom_bar(mapping = aes(x= Chromosome, fill = Species))
print(SNPs_per_chrom3)


SNP_density <- (ggplot(num_snp_fang, aes(x= as.numeric(Position)))) + geom_density(aes(fill = Chromosome)) + facet_wrap(~ Chromosome) +
                  labs(x = "Position", y = "Density") + ggtitle("SNP Density")
SNP_density2 <- (ggplot(genobound, aes(x= as.numeric(Position)))) + geom_density(aes(fill = Chromosome)) + facet_wrap(~ Chromosome) +
  labs(x = "Position", y = "Density") + ggtitle("SNP Density")
print(SNP_density)

add_column <- num_snp_fang
add_column$Heterozygotes <- "Heterozygotes"
add_column$Heterozygotes[add_column$NT == "?/?"] <- "Missing"
add_column$Heterozygotes[add_column$NT %in% c("A/A", "T/T", "C/C", "G/G")] <- "Homozygous"
summary(add_column)

hetero <- ggplot(add_column, aes(x = Sample_ID, fill = Heterozygotes))+ geom_bar(position = "fill")+ labs(x="Sample_ID", y="Proportion")+ggtitle("Heterozyosity in Miaze and Teo")
print(hetero)


Heterozygosity_of_All <- (ggplot(add_column, aes(x = Sample_ID, fill = Heterozygotes)) +
                            geom_bar(position = "fill") + labs(x = "Sample_ID", y = "Proportion") +
                            ggtitle("Heterozygosity in Maize and Teosinte") +
                            theme(plot.title = element_text(hjust = 0.5)))
print(Heterozygosity_of_All)
Heterozygosity_Across_Groups <- (ggplot(add_column, aes(x = Group, fill = Heterozygotes)) +
                                   geom_bar(position = "fill") +
                                   labs(x = "Group", y = "Proportion") + theme(axis.text.x = element_text(angle = 90)) +
                                   ggtitle("Heterozygosity for all Groups") + theme(plot.title = element_text(hjust = 0.5)))
print(Heterozygosity_Across_Groups)

added_column2 <- filter(add_column,Group == "ZMMIL" | Group == "ZMMLR" | Group == "ZMMMR")
My_Plot_Heterozygosity_In_Maize <- (ggplot(added_column2, aes(x = Group, fill = Heterozygotes)) +
                                      geom_bar(position = "fill") +
                                      labs(x = "Maize Group", y = "Proportion") +
                                      ggtitle("Heterozygosity in Just Maize") +
                                      theme(plot.title = element_text(hjust = 0.5)))
print(My_Plot_Heterozygosity_In_Maize)

Names <- colnames(geno)[-c(1:3)]
genotypes.melt <- melt(geno, measure.vars = Names)
colnames(genotypes.melt)[c(3,4,5)] <- c("Group","SNP_ID", "Allele") #Melting data to make it easier to work with
genotypes.melt$Ho <- (genotypes.melt$Allele =="A/A" | genotypes.melt$Allele =="C/C" | genotypes.melt$Allele =="G/G" | genotypes.melt$Allele =="T/T")
sortedmelt <- arrange(genotypes.melt, Sample_ID, Group)
summarize.ID <- ddply(sortedmelt, c("Sample_ID"), summarise, total_ho = sum (Ho, na.rm=TRUE), total_het = sum (!Ho, na.rm=TRUE), missing = sum(is.na(Ho)))#Missing data and heterozygous ratio parameters
summarize.melt <- melt(summarize.ID, measure.vars = c("total_ho", "total_het", "missing"))
ggplot(summarize.melt,aes(x = Sample_ID, y = value, fill=variable)) + geom_bar(stat = "identity", position = "stack")
attributes(summarize.melt)
