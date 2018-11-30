
#install pathview
source("http://bioconductor.org/biocLite.R")
install.packages("dplyr")
install.packages("tidyverse")
install.packages("cowplot")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pathview", version = "3.8")


#Load libraries
library(tidyverse)
library(cowplot)
library(pathview)
library(dplyr)



#Read in Data
KAAS_dat <- read_tsv(file="allMAGs_ko.cleaned.txt", col_names = FALSE)%>%
  dplyr::rename(Sequence_ID=X1, ko=X2)
checkM_dat <- read_tsv(file="SI_200m_checkM_stdout.tsv", col_names = TRUE)
bac120 <- read_tsv(file = "gtdbtk.bac120.classification_pplacer.tsv", col_names = FALSE)
ar122 <- read_tsv(file = "gtdbtk.ar122.classification_pplacer.tsv", col_names = FALSE)
chem_dat <- read_csv(file="Saanich_TimeSeries_Chemical.csv", col_names = TRUE)
allMAGs_dat <- read_tsv(file = "allMAGs.tsv", col_names = TRUE)
merged_MAGs_dat <- read_tsv(file= "merged_MAGs.tsv", col_names = TRUE)
KO_log_test <- read_csv(file="KO_log.csv", col_names = FALSE) %>% 
  separate(X1, into = c("1", "2", "3", "4", "5", "6")) %>% 
  dplyr::select("1", "3", "4") %>%
  dplyr::rename(metabolism = "1", ko = "3", gene = "4")
Nitrogen_KO_log <- read_csv("Nitrogen_KO_log.csv", col_names = FALSE) %>%
  dplyr::rename(ko = X1, gene = X2)
Sulfur_KO_log <- read_csv(file = "sulfur_ko_log.csv", col_names = FALSE) %>%
  dplyr::rename(ko = X1, gene = X2)
prokka_mag_log <- read_csv(file = "Prokka_MAG_log.csv", col_names = FALSE) %>%
  separate(X2, into = c("a", "b", "c", "d"), sep = "/") %>%
  select(mag = X1, Bin = c)

All_data_final <- left_join(All_data_final, KO_log_test)

View(KO_log_test)
#Transcriptional RPKM Data
SI042_rpkm_dat <- read_csv(file="SI042_allMAGs_RPKM.csv", col_names = FALSE) %>%
  dplyr::rename(Sequence_ID=X1, SI042_rpkm=X2)
SI048_rpkm_dat <- read_csv(file = "SI048_allMAGs_RPKM.csv", col_names = FALSE) %>%
  dplyr::rename(Sequence_ID=X1, SI048_rpkm=X2)
SI072_rpkm_dat <- read_csv(file = "SI072_allMAGs_RPKM.csv", col_names = FALSE) %>%
  dplyr::rename(Sequence_ID=X1, SI072_rpkm=X2)
SI073_rpkm_dat <- read_csv(file = "SI073_allMAGs_RPKM.csv", col_names = FALSE) %>%
  dplyr::rename(Sequence_ID=X1, SI073_rpkm=X2)
SI074_rpkm_dat <- read_csv(file = "SI074_allMAGs_RPKM.csv", col_names = FALSE) %>%
  dplyr::rename(Sequence_ID=X1, SI074_rpkm=X2)
SI075_rpkm_dat <- read_csv(file = "SI075_allMAGs_RPKM.csv", col_names = FALSE) %>%
  dplyr::rename(Sequence_ID=X1, SI075_rpkm=X2)

#Merge all cruises' RPKM Data
SI_ALL <- Reduce(function(x, y) merge(x,y, all =TRUE), list(SI042_rpkm_dat, SI048_rpkm_dat, SI072_rpkm_dat, SI073_rpkm_dat, SI074_rpkm_dat, SI075_rpkm_dat))

#Calculate the total RPKM across all cruises
rpkm_total <- mutate(SI_ALL, rpkm_sum = (apply(SI_ALL[2:7], 1, sum))) %>%
  select(Sequence_ID, rpkm_sum)


#Updated pathview test
all_ko_rpkm <- left_join(KAAS_dat, rpkm_total, by="Sequence_ID") %>%
  separate(Sequence_ID, into=c("mag", "orf")) %>%
  left_join(prokka_mag_log, by="mag") %>%
  left_join(Cleaned_dat, by="Bin")

crenarch_t_rpkm <- all_ko_rpkm %>%
  filter(Phylum == "p__Crenarchaeota") %>%
  group_by(mag, ko) %>%
  summarise(total = sum(rpkm_sum)) %>%
  spread(key = mag, value = total)

Crenarc_pv_mat <- dplyr::select(crenarch_t_rpkm, -ko)
rownames(Crenarc_pv_mat) <- crenarch_t_rpkm$ko

View(all_ko_rpkm)

write.csv(Crenarc_pv_mat, file = "Crenarc_pv_mat.csv")

#This does not work, error downloading xml and png files. Upliaded .csv file to online pathview portal instead
pv.out <- pathview(gene.data = pv_mat,
                   species = "ko",
                   pathway.id = "00910",
                   kegg.dir = "~/Desktop/MICB405_figures/")


#Add on the ko data
SI_ALL_ko <- merge(SI_ALL, KAAS_dat, by="Sequence_ID") %>% 
  separate(Sequence_ID, into=c("mag", "orf_id"))

View(KAAS_dat)

#Add gene data for each ko, but only for genes in the nitrogen AND sulfur pathways
SI_ALL_NS_ko <- merge(SI_ALL_ko, KO_log, by="ko") %>%
  dplyr::select(mag, ko, gene, SI042_rpkm, SI048_rpkm, SI072_rpkm, SI073_rpkm, SI074_rpkm, SI075_rpkm)

SI_rpkm_ko <- mutate(SI_ALL_NS_ko, total = (apply(SI_ALL_NS_ko[4:9], 1, sum)))

All_data_test <- merge(SI_rpkm_ko, prokka_mag_log, by="mag")
All_data_final <- merge(All_data_test, Cleaned_dat, by="Bin")

write.table(All_data_final, "All_data_final.csv", sep = ",")

####RPKM bubbleplot

n_genes <- All_data_final %>%
  filter(metabolism == "Nitrogen")
  
s_genes <- All_data_final %>%
  filter(metabolism == "Sulfur")
View(n_genes)

rpkm_plot <- ggplot() +
  geom_point(data=n_genes, aes(x = gene, y = mag, size = total, color = Phylum))+
  geom_point(data=s_genes, aes(x = gene, y = mag, size = total, color = Phylum))+
  labs(size = "Total RPKM") +
  scale_color_hue(labels = c("Bacteroidota", "Chloroflexota", "Crenarchaeota", "Elusimicrobiota", "Gemmatimonadota",
                                   "Marinisomatota", "Nanoarchaeota", "Patescibacteria", "Proteobactera", "SAR324")) +
  xlab("Metabolic Genes") +
  ylab("Metagenome Assembled Genomes") +
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=90, hjust=1, size = 8),
        axis.text.y=element_blank(),
        legend.position = "bottom") +
  facet_grid(~metabolism, scales = "free_x", space = "free_x")
  
rpkm_plot

####Sulfur Cycle Genes' Pie Charts

cysCDHN <- All_data_final %>%
  filter(gene == "cysC" | gene == "cysD" | gene == "cysH" | gene == "cysN", total < 100)
View(cysCDHN)

cysCDHN_pie_plot <- ggplot() +
  geom_bar(data = cysCDHN, stat = "identity", aes(x="", y=total, fill = Phylum)) +
  coord_polar("y", start = 0) +
  theme(axis.text.x=element_blank(),
        text = element_text(size=30)) +
  scale_fill_manual(labels = c("Bacteroidota", "Chloroflexota", "Crenarchaeota", "Patescibacteria", "Proteobacteria"),
                    values = c("#F8766D", "#DE8C00", "#B79F00", "#00B4F0", "#F564E3")) +
  xlab(element_blank()) +
  ylab(element_blank())
cysCDHN_pie_plot

##

cysJI_sir <- All_data_final %>%
  filter(str_detect(gene, "sir") | gene == "cysJ" | gene == "cysI")
View(cysJI_sir)

cysJI_sir_pie_plot <- ggplot() +
  geom_bar(data = cysJI_sir, stat = "identity", aes(x="", y=total, fill = Phylum)) +
  coord_polar("y", start = 0) +
  theme(axis.text.x=element_blank(),
        text = element_text(size=30)) +
  scale_fill_manual(labels = c("Bacteroidota", "Proteobacteria", "SAR324"),
                    values = c("#F8766D", "#F564E3", "#FF64B0" )) +
  xlab(element_blank()) +
  ylab(element_blank())
cysJI_sir_pie_plot         

##


apr <- All_data_final %>%
  filter(str_detect(gene, "apr") | str_detect(gene, "sat") | str_detect(gene, "dsr"), total < 100)
View(apr)

apr_pie_plot <- ggplot() +
  geom_bar(data = apr, stat = "identity", aes(x="", y=total, fill = Phylum)) +
  coord_polar("y", start = 0) +
  theme(axis.text.x=element_blank(),
        text = element_text(size=30)) +
  scale_fill_manual(labels = c("Bacteroidota", "Chloroflexota", "Crenarchaeota", "Marinisomatota", "Proteobacteria", "SAR324"),
                    values = c("#F8766D", "#DE8C00", "#B79F00", "#00C08B" , "#F564E3", "#FF64B0")) +
  xlab(element_blank()) +
  ylab(element_blank())
apr_pie_plot

##

soe <- All_data_final %>%
  filter(str_detect(gene, "soe"), total < 100)
View(soe)

##

sqr <- All_data_final %>%
  filter(str_detect(gene, "sqr"), total <100)
View(sqr)

sqr_pie_plot <- ggplot() +
  geom_bar(data = sqr, stat = "identity", aes(x="", y=total, fill = Phylum)) +
  coord_polar("y", start = 0) +
  theme(axis.text.x=element_blank(),
        text = element_text(size=30)) +
  scale_fill_manual(labels = c("Bacteroidota", "Chloroflexota", "Proteobacteria", "SAR324"),
                    values = c("#F8766D", "#DE8C00", "#F564E3", "#FF64B0" )) +
  xlab(element_blank()) +
  ylab(element_blank())
sqr_pie_plot 

sox <- All_data_final %>%
  filter(str_detect(gene, "sox"))
View(sox)


####Nitrogen Cycle Genes' Pie Charts

View(All_data_final)

#Filter the data for each gene
nap_nar <- All_data_final %>%
  filter(gene == "napA" | gene == "narG", total < 100)
View(nap_nar)

nar <- All_data_final %>%
  filter(str_detect(gene, "nar"), total < 100)
View(nar)

nrfA <- All_data_final %>%
  filter(str_detect(gene, "nrfA"))
View(nrfA)

hao_amo <- All_data_final %>%
  filter(str_detect(gene, "hao") | str_detect(gene, "amo"))
View(hao_amo)

gdh_gln <- All_data_final %>%
  filter(str_detect(gene, "gdh") | str_detect(gene, "gln"), total < 100)
View(gdh_gln)

nir <- All_data_final %>%
  filter(str_detect(gene, "nir"), total < 100)
View(nir)

nirs <- All_data_final %>%
  filter(str_detect(gene, "nirS"))

View(nirs)
nor <- All_data_final %>%
  filter(str_detect(gene, "nor"))
View(nor)

nos <- All_data_final %>%
  filter(str_detect(gene, "nos"))
View(nos)  

library(scales)
show_col(hue_pal()(12))


#Generate the nitrogen plots
nap_nar_pie_plot <- ggplot() +
  geom_bar(data = nap_nar, stat = "identity", aes(x="", y=total, fill = Phylum)) +
  coord_polar("y", start = 0) +
  theme(axis.text.x=element_blank()) +
  scale_fill_manual(labels = c("Chloroflexota", "Marinisomatota", "Proteobacteria"),
                     values = c("#DE8C00", "#00BFC4", "#F564E3")) +
  xlab(element_blank()) +
  ylab(element_blank())
 
nar_pie_plot <- ggplot() +
  geom_bar(data = nar, stat = "identity", aes(x="", y=total, fill = Phylum)) +
  coord_polar("y", start = 0) +
  theme(axis.text.x=element_blank()) +
  scale_fill_manual(labels = c("Chloroflexota", "Marinisomatota", "Proteobacteria"),
                    values = c("#DE8C00", "#00BFC4", "#F564E3")) +
  xlab(element_blank()) +
  ylab(element_blank())

nir_pie_plot <- ggplot() +
  geom_bar(data = nir, stat = "identity", aes(x="", y=total, fill = Phylum)) +
  coord_polar("y", start = 0) +
  theme(axis.text.x=element_blank()) +
  scale_fill_manual(labels = c("Bacteroidota", "Crenarchaeota", "Proteobacteria"),
                    values = c("#F8766D", "#B79F00", "#F564E3")) +
  xlab(element_blank()) +
  ylab(element_blank())

nor_pie_plot <- ggplot() +
  geom_bar(data = nor, stat = "identity", aes(x="", y=total, fill = Phylum)) +
  coord_polar("y", start = 0) +
  theme(axis.text.x=element_blank()) +
  scale_fill_manual(labels = c("Proteobacteria"),
                    values = c("#F564E3")) +
  xlab(element_blank()) +
  ylab(element_blank())

nos_pie_plot <- ggplot() +
  geom_bar(data = nos, stat = "identity", aes(x="", y=total, fill = Phylum)) +
  coord_polar("y", start = 0) +
  theme(axis.text.x=element_blank()) +
  xlab(element_blank()) +
  ylab(element_blank())

nirs_pie_plot <- ggplot() +
  geom_bar(data = nirs, stat = "identity", aes(x="", y=total, fill = Phylum)) +
  coord_polar("y", start = 0) +
  theme(axis.text.x=element_blank()) +
  xlab(element_blank()) +
  ylab(element_blank())


nrfA_pie_plot <- ggplot() +
  geom_bar(data = nrfA, stat = "identity", aes(x="", y=total, fill = Phylum)) +
  coord_polar("y", start = 0) +
  theme(axis.text.x=element_blank()) +
  scale_fill_manual(labels = c("Not Found"),
                    values = c("#F8766D")) +
  xlab(element_blank()) +
  ylab(element_blank())

gdh_gln_pie_plot <- ggplot() +
  geom_bar(data = gdh_gln, stat = "identity", aes(x="", y=total, fill = Phylum)) +
  coord_polar("y", start = 0) +
  theme(axis.text.x=element_blank(),
        text = element_text(size=30)) +
  scale_fill_manual(labels = c("Bacteroidota", "Chloroflexota", "Crenarchaeota", "Marinisomatota", "Nanoarchaeota", "Patescibacteria", "Proteobacteria"),
                    values = c("#F8766D", "#DE8C00", "#B79F00", "#00C08B", "#00BFC4", "#00B4F0", "#F564E3")) +
  xlab(element_blank()) +
  ylab(element_blank())

gdh_gln_pie_plot
nrfA_pie_plot
nar_pie_plot
nirs_pie_plot
nap_nar_pie_plot
nir_pie_plot
nor_pie_plot
nos_pie_plot

#Contamination vs Completion Figure

merged_MAGs_dat <- dplyr::rename(merged_MAGs_dat, Bin=`project2/prokka_output/SaanichInlet_200m.107/MAG107.tsv`)
  
merged_MAGs_dat <- merge(merged_MAGs_dat, KAAS_dat, by="locus_tag", all=TRUE)



###bin quality figure

#first select and filter the relevant data
bin_dat <- checkM_dat %>%
  dplyr::select(`Bin Id`, `Marker lineage`, `Completeness`, `Contamination`) %>%
  filter(Contamination <30) %>%
  dplyr::rename(Bin = `Bin Id`)
  
View(bin_dat)

#Rename column Names for classification data
bac120 <- dplyr::rename(bac120, Bin = `X1`)
bac120 <- dplyr::rename(bac120, Classification = X2)

ar122 <- dplyr::rename(ar122, Bin = `X1`)
ar122 <- dplyr::rename(ar122, Classification = X2)

#merge classification data to bin data and separate each clade into new columns

merged_class <- rbind(bac120, ar122)
merged_dat <- merge(bin_dat, merged_class, by="Bin", all=TRUE)
Cleaned_dat <- separate(data = merged_dat, col = Classification, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"), sep = "\\;")

View(Cleaned_dat)


#Generate a scatterplot to visualize Contaminate vs Completeness, with Clades colored (size to be done...)
bin_fig <- Cleaned_dat %>%
  ggplot() +
  geom_point(aes(x=Completeness, y=Contamination, color=Phylum))+
  xlab("Completeness (%)") +
  ylab("Contamination (%)") +
  theme_grey(base_size = 12) +
  theme(legend.text = element_text(size = 8), legend.position = "bottom")
  
bin_fig

