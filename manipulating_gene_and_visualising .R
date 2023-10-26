#script manipulate the gene experssion data
# setwd("~/Downloads/bio R project")

install.packages("tidyverse")
#load libraries
library(dplyr)
library(tidyverse)
library(GEOquery)
library(ggplot2)

#read in the data
 dat <- read.csv("GSE183947_fpkm.csv")
dim(dat)

#get metadata
gse <- getGEO(GEO = 'GSE183947', GSEMatrix = TRUE)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10000)

gse 

metadata <- pData(phenoData(gse[[1]]))

head(metadata)
metadata.subset <- select(metadata, c(1,10,11,17))

metadata.modefied <- metadata %>%
  select(1,10,11,17) %>%
  rename(tissue = characteristics_ch1) %>%
  rename(metastasis = characteristics_ch1.1) %>%
  mutate(tissue = gsub("tissue: ", "",tissue))%>%
  mutate(metastasis = gsub("metastasis: ", "",metastasis))
  
head(dat)

#reshaping data

data.long <- dat %>%
  rename(gene = X)%>%
  gather(key = 'samples', value = 'FPKM', -gene) 

#join dataframes = dat.long + metadata.modified

data.long <- data.long %>%
  left_join(., metadata.modefied, by = c("samples" = "description"))

#explore the data

data.long %>%
  filter(gene == "BRCA1" | gene == "BRCA2") %>%
  group_by(gene, tissue) %>%
  summarise(mean_FPKM = mean(FPKM))%>%
  head()

#visualise the gene

# 1 barplot
data.long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = samples, y = FPKM, fill = tissue)) +
  geom_col()

# 2 density 

data.long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = FPKM, fill = tissue)) +
  geom_density(alpha = .5)

# 3 boxplot
data.long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = metastasis, y = FPKM)) +
  geom_boxplot() 


# scatter pot
data.long %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
  spread(key = gene, value = FPKM) %>%
  ggplot(., aes(x= BRCA1, y=BRCA2, color = tissue))+
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE)


#heatmap

genes.of.interest <- c('BRCA1', 'BRCA2', 'TP53', 'ALK', 'MYCN')

data.long%>%
  filter(gene %in% genes.of.interest)%>%
  ggplot(., aes(x=samples, y=gene, fill = FPKM))+
  geom_tile() +
  scale_fill_gradient(low = 'white', high = 'red')














