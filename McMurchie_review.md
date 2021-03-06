# Elizabeth McMurchie suggestions for Nic Lizer's R project
### Specific errors
- The code chunk starting with line 85 has a serious, code-breaking error where the object `snp_maizegeno.select` is not found. You must have missed a step where you specified how this object was defined. This would explain why none of your 40 files made it into the final directory that I pulled from Github. Assign the right values to this object, and things should work out!
- Note: I ran into the same issue with line 182, so I wasn't able to reproduce your graph.
- I had an issue with the second graph, as well: `Error in add_column$value : object of type 'closure' is not subsettable`. I'm not entirely sure what this means, but since the issue appears to be with creating columns for homozygosity, heterozygosity, etc. you might try something like I did:
- ```geno5 <- geno %>%
  mutate(homozygous = apply(geno, 1, function(x) length(which(x == ("A/A") | x == ("T/T") | x == ("C/C") | x == ("G/G"))))) %>%
  mutate(heterozygous = apply(geno, 1, function(x) length(which(x == ("A/C") | x == ("A/T") | x == ("A/G") | x == ("C/A") | x == ("C/T") | x == ("C/G") | x == ("G/A") | x == ("G/C") | x == ("G/T") | x == ("T/A") | x == ("T/C") | x == ("T/G"))))) %>%
  mutate(missing = apply(geno, 1, function(x) length(which(x == ("?/?"))))) %>%
  mutate(total =  apply(geno, 1, function(x) length(which(x == ("?/?") | x == ("A/A") | x == ("T/T") | x == ("C/C") | x == ("G/G") | x == ("A/C") | x == ("A/T") | x == ("A/G") | x == ("C/A") | x == ("C/T") | x == ("C/G") | x == ("G/A") | x == ("G/C") | x == ("G/T") | x == ("T/A") | x == ("T/C") | x == ("T/G")))))```
You will then want to filter by group and select by Sample_id, Group, and your homozygosity, heterozygosity, and missing data columns. Then you'll need to pivot the homozygosity, heterozygosity, and missing data columns longer, so you have a "zygosity" (or similar - this may not be the best name) column and "counts" column.

### General suggestions
- Once you get your `snp+maizegeno.select` object assigned correctly, you should be able to write your files as part of a loop using something like:
```
 for (i in 1:10) {
  temp_df <- snp_maizegeno.select %>%
    filter(Chromosome == i) 
write.table(temp_df, file = paste("./directory_name/maize_chr",i,".tsv", sep = ""), row.names = FALSE, col.names = TRUE)
```
In this example I made tab-delimited files rather than comma-delimited files, but I don't think the distinction ultimately matters. 
- You generally don't want to write out your own filepath when saving things (as you did in line 137 onward). This isn't replicable for other users! Instead, you can use something like `write.csv(maize_chrom_inc8,"./maize_chrom_inc8")` which will put them into the main directory associated with this R project. It may also be helpful to have R print your files into subdirectories using `ggsave("./directory_name/plot_name.png", plot = plot_name etc.)` for the sake of organization. That's also what the `"//director_name..."` part of the code snippet above does.
- Say "knit" (the ball of thread and knitting needles on R studio above your script window) to make an HTML document that can be pushed online and easily read. 
- You spelled "teosinte" incorrectly for some of your file names, which may be confusing. This error starts on line 124.

### General comments
- I can't really comment on the graphs because I wasn't able to run them and they weren't pushed to your Github. Let me know if you get these working and I can take a look!