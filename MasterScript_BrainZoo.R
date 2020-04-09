library(data.table)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(bsseq)
library(R.utils)
library(reshape2)
library(ggseqlogo)
library(cowplot)

######################################################################
## Preprocessing of data
######################################################################

# WGBS fastq files mapped to reference genomes using BSseeker2 with bowtie2 aligner / fastp trimmer / sambamba markdup
# Aligned bams are then processed with CGmapTools
# ATCGmap is parsed using ATCGmap_GenotypeCorrector.sh in this same repository.


######################################################################
## Reading CGmap into R and export RDS with bsseq objects (CG)
######################################################################

options(scipen=999)

# Set command line arguments
args <- commandArgs(TRUE)



CGmap <- args[1]
name <- args[2]
fai <- args[3]

fai <- fread(file = fai)[,c(1,2)] %>% .[!.$V1 %in% c("chrL","chrP"),]
sLengths <- fai$V2
names(sLengths) <- fai$V1


read_CGmap_into_CG_lambda_files <- function(CGmap, name, sLengths){
  
  # Read the file
  dat <- fread(input = CGmap, sep = "\t", select = c(1,2,3,4,5,7,8),
               col.names = c("chr", "base", "position", "context","dinucleotide",
                             "C_reads", "CT_reads"))
  # Subset to CG context only
  message(paste0("processing CG..."))
  datCG <- dat[dat$context == "CG" & !(dat$chr %in% c("chrL","chrP","MT","chrMT")), ]
  
  datCG$strand <- ifelse(test = datCG$base == "G",
                         yes = "-",
                         no = "+")
  
  bs_obj <- BSseq(gr = GRanges(seqnames = datCG$chr,
                               ranges = IRanges(start = as.numeric(datCG$position),
                                                end = as.numeric(datCG$position)), strand = datCG$strand , seqlengths = sLengths), 
                  sampleNames = name, M = as.matrix(datCG$C_reads), Cov = as.matrix(datCG$CT_reads),
                  rmZeroCov = FALSE)
  
  bs_obj_collapsed <- strandCollapse(bs_obj)
  rm(datCG)
  
  saveRDS(object = bs_obj_collapsed, file = paste0(name,".CG_bsseq.rds"))
  
  ###########
  message(paste0("processing Lambda..."))
  # Subset chrL
  dat_Lambda <- dat[dat$chr == "chrL", ]
  
  #add_sp_name
  dat_Lambda$species <- name
  
  #save lambda genome
  saveRDS(object = dat_Lambda, file = paste0(name,".lambda.rds"))
  
  gc()
  
  # return the aggregated data object
  message(paste0("Objects created for ",name))
}

read_CGmap_into_CG_lambda_files(CGmap = CGmap, name = name, sLengths = sLengths)

######################################################################
### Reading lambda genomes and get non-conversion rates
######################################################################


Lambda_to_nonConversion <- function(bsobj){
  a <- readRDS(bsobj)
  name <- str_split(bsobj, pattern = "\\.", simplify = T)[1] 
  b <- a %>% data.frame() %>% group_by(.,dinucleotide) %>% summarise(sum(C_reads),sum(CT_reads),n())
  colnames(b) <- c("context","mC","C","positions")
  b <- mutate(b, mC_C = 100*mC/C)
  message("on")
  c <- data.frame(context = "CH", mC = sum(b$mC), C = sum(b$C), positions = sum(b$positions)) %>% mutate(mC_C = 100*mC/C)
  b <- rbind(b,c)
  message("peta")
  b$species <- name
  b$variable <- paste0("Lambda_m",b$context)
  b <- b[b$context != "--",]
  return(b)
}

list_of_files <- list.files(pattern = ".lambda.rds$", full.names = TRUE) %>% str_replace(pattern = "\\./", replacement = "")

non_conversion_rates <- lapply(list_of_files, FUN = Lambda_to_nonConversion) %>% do.call(rbind,.) %>% data.frame()

write.table(x = non_conversion_rates, file = "non_conversion_rates.tsv", quote = F, sep = "\t")


######################################################################
# Read CGmap and establish CH methylation stats globally
######################################################################

read_CGmap_into_mCH_levels <- function(CGmap){
  
  # Read the file
  dat <- fread(input = CGmap, sep = "\t", select = c(1,2,3,4,5,7,8),
               col.names = c("chr", "base", "position", "context","dinucleotide",
                             "C_reads", "CT_reads"))
  
  ###########
  message(paste0("processing CH..."))
  # Subset CH
  dat <- dat[!grepl(dat$context,pattern = "CG") & !(dat$chr %in% c("chrL","chrP","chrMT","MT")), ]
  
  total_CHcov <- sum(dat$CT_reads)
  total_mCH <- sum(dat$C_reads)
  mCA_cov <- dat[dat$dinucleotide == "CA",]
  mCA_M <- sum(mCA_cov$C_reads)
  mCA_cov <- sum(mCA_cov$CT_reads)
  mCT_cov <- dat[dat$dinucleotide == "CT",]
  mCT_M <- sum(mCT_cov$C_reads)
  mCT_cov <- sum(mCT_cov$CT_reads)
  mCC_cov <- dat[dat$dinucleotide == "CC",]
  mCC_M <- sum(mCC_cov$C_reads)
  mCC_cov <- sum(mCC_cov$CT_reads)
  df <- data.frame(file = CGmap, total_CHcov = total_CHcov, total_mCH = total_mCH,
                   mCA_cov = mCA_cov, mCA_M = mCA_M,
                   mCT_cov = mCT_cov, mCT_M = mCT_M,
                   mCC_cov = mCC_cov, mCC_M = mCC_M) %>% mutate(mCA = mCA_M/mCA_cov,mCT = mCT_M/mCT_cov,mCC = mCC_M/mCC_cov)
  
  
  rm(dat)
  gc()
  
  # return the aggregated data object
  return(df)
}

read_CGmap_into_mCH_levels(CGmap = CGmap)

list_of_files <- list.files(pattern = "corrected.gz$", full.names = TRUE) %>% str_replace(pattern = "\\./", replacement = "")

mCH_table <- mclapply(list_of_files,FUN = read_CGmap_into_mCH_levels, mc.cores = 3) %>% do.call(rbind,.) %>% data.frame(stringsAsFactors = F)

write.table(x = mCH_table, file = "mCH_stats.tsv", quote = F, sep = "\t")

######################################################################
# Compare CH levels to lambda non-conversion levels
######################################################################

#mCH_table <-fread("mCH_stats.tsv") %>% data.frame() %>% dplyr::select(-V1)
#non_conversion <- fread("non_conversion_rates.tsv") %>% data.frame() %>% dplyr::select(-V1)

mCH_table$species <- factor(c("Shark","Gallus","Lamprey","Mouse","Opossum","Platypus","Octopus","Xenopus","Zebrafish","Honeybee","Human","Parus","Zebrafish_forebrain","Lancelet_neural"),
                           levels = rev(c("Human", "Mouse","Opossum","Platypus","Parus","Gallus","Xenopus","Zebrafish_forebrain","Shark","Lamprey","Lancelet_neural","Octopus","Honeybee")))


mCH_to_plot <- mCH_table %>% dplyr::select(species, mCA, mCT, mCC) %>% melt() %>% mutate(value = value*100, genome = "Nuclear")

non_conversion_toPlot <- non_conversion %>% mutate(variable = str_replace(variable, pattern = "Lambda_", replacement = "")) %>% 
  dplyr::rename(value = mC_C) %>% dplyr::select(species, variable, value) %>% filter(!variable %in% c("mCG","mCH")) %>%
  mutate(genome = "Lambda")

mCH_context_to_plot <- rbind(mCH_to_plot,non_conversion_toPlot )
mCH_context_to_plot$genome <- factor(mCH_context_to_plot$genome, levels = rev(c("Nuclear","Lambda")))

gg_barplot_CH_nonconversion <- ggplot(mCH_context_to_plot, aes(x = species, y = value, fill = genome)) + geom_bar(stat="identity", width=.5, position = "dodge") + 
  facet_grid(. ~ variable) + ylab("mC/C methylation %") + scale_fill_brewer(palette="Paired") + coord_flip()

ggsave(filename = "Barplot_CHcontext_nonconversion.pdf", plot = gg_barplot_CH_nonconversion)


######################################################################
## Load CG rds and get the total mCG and plot
######################################################################

mCG_on_CG <- function(bsobj){
  bs <- readRDS(bsobj)
  df <- data.frame(mCG = sum(getCoverage(bs, type = "M")[,1])/sum(getCoverage(bs, type = "Cov")[,1]),
                   species = str_split(bsobj, pattern = "\\.", simplify = T)[1],
                   median_cov = median(getCoverage(bs, type = "Cov")[,1]),
                   CG_cov = sum(getCoverage(bs, type = "Cov")[,1]), CG_M = sum(getCoverage(bs, type = "M")[,1]))
  return(df)
}

list_of_files <- list.files(pattern = ".CG_bsseq.rds$", full.names = TRUE) %>% str_replace(pattern = "\\./", replacement = "")

mCG_total_df <- mclapply(list_of_files,FUN = mCG_on_CG, mc.cores = 4) %>% do.call(rbind,.) %>% data.frame(stringsAsFactors = F)

write.table(x = mCG_total_df, file = "Multispecies_mCG_table.tsv", quote = F, sep = "\t")

#Read mCG table and plot

mCG_total_df <- fread("Multispecies_mCG_table.tsv") %>% data.frame() %>% mutate(mCG_perc = mCG *100)
mCG_total_df$species <- factor(mCG_total_df$species, levels = c(levels = rev(c("Human","Mouse","Opossum","Platypus","Parus","Gallus","Xenopus","Zebrafish_forebrain",
                                                                               "Shark","Lamprey",
                                                                               "Ciona","Lancelet_neural","Octopus","Honeybee","Nematostella"))))


gg_mCG <- ggplot(mCG_total_df, aes(x = species, y = mCG_perc)) + geom_bar(stat= "identity") + theme_bw() + ylab("mCG/CG Global Methylation %") +
  theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, hjust = 1)) + ylim(c(0,100)) + coord_flip()

ggsave(filename = "mCG_barplot.pdf", plot = gg_mCG)


######################################################################
## Load CG rds and get the OZ global mC proportion plots
######################################################################


mCG_to_Oz_proportions <- function(bsobj){
  bs <- readRDS(bsobj)
  
  #filter for 10x coverage
  ten_x_covered_CGs <- getMeth(bs, type = "raw")[getCoverage(bs, type = "Cov") >= 10, 1]
  rm(bs)
  gc()
  
  
  #count positions
  none <- sum(ten_x_covered_CGs == 0) / length(ten_x_covered_CGs)
  low <- sum(ten_x_covered_CGs >0 & ten_x_covered_CGs <=0.2) / length(ten_x_covered_CGs)
  inter <- sum(ten_x_covered_CGs >0.2 & ten_x_covered_CGs <=0.8) / length(ten_x_covered_CGs)
  high <- sum(ten_x_covered_CGs >0.8) / length(ten_x_covered_CGs)
  
  df <- data.frame(type = c("none","low","inter","high"),
                   species = str_split(bsobj, pattern = "\\.", simplify = T)[1],
                   proportion = c(none, low, inter, high)) %>%
    mutate(perc = round(proportion*100, digits = 2))
  return(df)
}

list_of_files <- list.files(pattern = ".CG_bsseq.rds$", full.names = TRUE) %>% str_replace(pattern = "\\./", replacement = "")

OzPlotProportions <- mclapply(list_of_files, FUN = mCG_to_Oz_proportions,mc.cores = 5) %>% do.call(rbind,.) %>% data.frame()

OzPlotProportions$type <- factor(OzPlotProportions$type, levels = c("none","low","inter","high"))

species <- c("Nematostella","Honeybee","Octopus",
             "Lancelet_neural","Ciona","Lamprey","Shark",
             "Zebrafish_forebrain","Xenopus","Gallus","Parus","Platypus",
             "Opossum","Mouse","Human")
OzPlotProportions <- OzPlotProportions[OzPlotProportions$species %in% species,]
OzPlotProportions$species <- factor(OzPlotProportions$species, levels = c("Nematostella","Honeybee","Octopus",
                                                                          "Lancelet_neural","Ciona","Lamprey","Shark",
                                                                          "Zebrafish_forebrain","Xenopus","Gallus","Parus","Platypus",
                                                                          "Opossum","Mouse","Human"))


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

gg_Oz <- ggplot(data = OzPlotProportions, aes(x = as.factor(species),y = perc, fill = type)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous() + 
  # X axis label
  xlab(label = "") + 
  # Y axis label
  ylab(label = "mCG/CG") +
  #theme with white background
  theme_bw() +
  #eliminates background, gridlines, and chart border
  theme( axis.text.x = element_text(angle = 90, hjust = 1, size = 15)
         ,plot.background = element_blank()
         ,panel.grid.major = element_blank()
         ,panel.grid.minor = element_blank()
         ,panel.border = element_blank()
  ) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black')) +
  # Edit legend
  scale_fill_manual(values = cbPalette, name="mCG/CG level",
                    labels=c("No (mCG/CG = 0) ",
                             "Low (mCG/CG >0 and <0.2)",
                             "Intermediate (mCG/CG >0.2 and <0.8)",
                             "High (mCG/CG > 0.8)"))


ggsave(filename = "Multispecies_OZplot.pdf", plot = gg_Oz)

######################################################################
# Get the motifs for top 10,000 mCH positions
######################################################################

# CGmap files parsed with bedtools and bash:
#mC=0.1                 # minimal mC/C level
#cov=10                 # minimal coverage required
#max_positions=10000    # total number of positions
#genome=$1              # genome fasta file
#CGmap=$2               # CGmap file
#genome_fai=$(echo $genome |awk '{print $0".fai"}')
#name=$(echo $CGmap|perl -pe "s/\.corrected.gz//" |perl -pe "s/\.CGmap//")
#zcat $CGmap |awk -v cove="$cov" -v mClev="$mC" '$4 ~ /CH/ && $8 >= cove && $6 >= mClev' |grep -vE "chrM|MT|chrL|chrP" | \
#bioawk -t '{if($2 == "C"){print $1,$3-1,$3,$4,$6,"+"}else if($2 == "G"){print $1,$3-1,$3,$4,$6,"-"}}' | \
#sort -k5,5nr |head -n $max_positions |sort -k1,1V -k2,2n -k3,3n > ${name}_${mC}_${cov}.bed;
#bedtools slop -i ${name}_${mC}_${cov}.bed -g $genome_fai -b 5 > ${name}_${mC}_${cov}.5flank.bed;
#bedtools getfasta -bed ${name}_${mC}_${cov}.5flank.bed -fi $genome -fo ${name}_${mC}_${cov}.5flank.fasta -s;
#seqkit fx2tab ${name}_${mC}_${cov}.5flank.fasta|cut -f 2 > ${name}_${mC}_${cov}.5flank.seqs


#seq Logo from Bedtools parsing
dir <- "~/TopCH/"

mCH_seqs <- list(scan(file = paste0(dir,"Bee_brain_lister_0.1_10.5flank.seqs"), what = "character"),
                 scan(file = paste0(dir,"RL125_Octopus_supraesophagealbrain_mC_merge_0.1_10.5flank.seqs"), what = "character"),
                 scan(file = paste0(dir,"MethC_Amphioxus_neural_tube_both_runs_0.1_10.5flank.seqs"), what = "character"),
                 scan(file = paste0(dir,"lamprey_brain_0.1_10.5flank.seqs"), what = "character"),
                 scan(file = paste0(dir,"eshark_0.1_10.5flank.seqs"), what = "character"),
                 scan(file = paste0(dir,"Zebrafish_forebrain_0.1_10.5flank.seqs"), what = "character"),
                 scan(file = paste0(dir,"xenopus_brain_0.1_10.5flank.seqs"), what = "character"),
                 scan(file = paste0(dir,"Gallus_gallus_5_brain_methylome_0.1_10.5flank.seqs"), what = "character"),
                 scan(file = paste0(dir,"Parus_0.1_10.5flank.seqs"), what = "character"),
                 scan(file = paste0(dir,"Platypus_brain_WGBS_0.1_10.5flank.seqs"), what = "character"),
                 scan(file = paste0(dir,"Opossum_brain_methylome_0.1_10.5flank.seqs"), what = "character"),
                 scan(file = paste0(dir,"MethylC-Seq_mm_fc_6wk_0.1_10.5flank.seqs"), what = "character"),
                 scan(file = paste0(dir,"MethylC-Seq_hs_mfg_25yr_0.1_10.5flank.seqs"), what = "character")
)

names(mCH_seqs) <- c("Honeybee", "Octopus", "Lancelet", "Lamprey", "Elephant shark","Zebrafish_forebrain",
                     "Zebrafish", "Xenopus", "Chicken", "Parus","Platypus","Opossum", "Mouse", "Human")

add_nums_to_names <- function(x){
  name <- names(mCH_seqs[x])
  num <- length(mCH_seqs[[x]])
  name <- paste(name, " ( n = ",num,")", sep = "")
  return(name)
}

names(mCH_seqs) <- lapply(1:length(mCH_seqs), add_nums_to_names) %>% unlist()

gg_mCH_logos <- ggseqlogo( mCH_seqs,  seq_type='dna',  ncol = 2,method = 'bits')

ggsave(filename = "mCH_logos.pdf", plot = gg_mCH_logos)

######################################################################
## Reading CGmap into R and export RDS with bsseq objects (CpA)
######################################################################

# Set command line arguments
args <- commandArgs(TRUE)
path <- "~/BrainZoo_RDS/"

CGmap <- args[1]
name <- args[2]
fai <- args[3]

fai <- fread(file = fai)[,c(1,2)] %>% .[!.$V1 %in% c("chrL","chrP"),]
sLengths <- fai$V2
names(sLengths) <- fai$V1


read_CGmap_into_CA_bsseq <- function(CGmap, name, sLengths){
  
  # Read the file
  dat <- fread(input = CGmap, sep = "\t", select = c(1,2,3,4,5,7,8),
               col.names = c("chr", "base", "position", "context","dinucleotide",
                             "C_reads", "CT_reads"))
  # Subset to CG context only
  message(paste0("processing CA..."))
  dat <- dat[dat$dinucleotide == "CA" & !grepl(dat$context,pattern = "CG") & !(dat$chr %in% c("chrL","chrP","MT","chrMT")), ]
  
  
  bs_obj <- BSseq(gr = GRanges(seqnames = dat$chr,
                               ranges = IRanges(start = as.numeric(dat$position),
                                                end = as.numeric(dat$position)), strand = dat$strand , seqlengths = sLengths), 
                  sampleNames = name, M = as.matrix(dat$C_reads), Cov = as.matrix(dat$CT_reads),
                  rmZeroCov = FALSE)
  
  saveRDS(object = bs_obj, file = paste0(path,name,".CA_bsseq.rds"))
  
  rm(dat,bs_obj)
  gc()
  
  #gc()
  
  # return the aggregated data object
  message(paste0("Objects created for ",name))
}

read_CGmap_into_CA_bsseq(CGmap = CGmap, name = name, sLengths = sLengths)


######################################################################
## Reading CpA bsseq, CpG bsseq, TPM file and gene annotation bed file
## and get mCH/mCG/TPM on genes
######################################################################

# Set command line arguments
args <- commandArgs(TRUE)

gene_bed <- args[1]
name <- args[2]
chr_append_logical <- args[3]

minimal_num_CpG = 30
minimal_mean_cov = 4

read_bed6_to_GRobject <- function(bedfile, chr_append = chr_append_logical){
  dat <- fread(bedfile)
  if (chr_append == T){
    dat$V1 <- paste0("chr",dat$V1)
  }
  #dat <- dat[dat$V1 %in% names(sLengths),]
  if (sum(grepl(unique(dat$V6), pattern = "C")) == 1){
    dat$V6 <- ifelse(dat$V6 == "C", yes = "-","+")
  }
  gr <- GRanges(seqnames = Rle(dat$V1),
                ranges = IRanges(start = dat$V2, end = dat$V3),
                strand = dat$V6, 
                gene_id = dat$V5, target_id = dat$V4)
  return(gr)
}

message(paste0("Reading ",name," genes..."))
gr <- read_bed6_to_GRobject(bedfile = gene_bed)

message(paste0("Reading ",name," TPMs..."))
tpm <- fread(paste0(dir,name,".kallisto.tsv")) %>% data.frame()
colnames(tpm) <- c("target_id","length","eff_length","est_counts","tpm")

message(paste0("Reading ",name," mCG..."))
bs_obj <- readRDS(paste0(name,".CG_bsseq.rds"))

coverage_positions <- getCoverage(bs_obj, regions = gr, type = "Cov", what = "perBase")
position_per_gene <- lapply(coverage_positions, FUN = length) %>% unlist()
mean_coverage_positions <- getCoverage(bs_obj, regions = gr, type = "Cov", what = "perRegionAverage")[,1]
M <- getCoverage(bs_obj, regions = gr, type = "M", what = "perRegionTotal")
C <- getCoverage(bs_obj, regions = gr, type = "Cov", what = "perRegionTotal")
mCG_on_regions <- (M/C)[,1]
rm(M,C)

message(paste0("Reading ",name," mCA..."))
bs_obj <- readRDS(paste0(name,".CA_bsseq.rds"))
M <- getCoverage(bs_obj, regions = gr, type = "M", what = "perRegionTotal")
C <- getCoverage(bs_obj, regions = gr, type = "Cov", what = "perRegionTotal")
mCA_on_regions <- (M/C)[,1]
rm(bs_obj, M, C)
gc()

df <- data.frame(target_id = gr$target_id, gene_length = width(gr), covered_Cs = position_per_gene, mean_coverage = mean_coverage_positions, 
                 mCG = mCG_on_regions, mCA = mCA_on_regions ) %>%
  left_join(.,tpm)

df <- df[df$covered_Cs > minimal_num_CpG & df$mean_coverage >= minimal_mean_cov & !is.na(df$tpm),]

#classify genes according to deciles of expression! 10 deciles and <1 TPM
message(paste0("Classifying..."))
df <- df %>% dplyr::select(-rank)
df$decile <- cut(df$tpm, breaks = c(0,quantile(df$tpm[df$tpm >= 1], seq(0, 1, length = 11), type = 5)),
                 labels = c("No expr.","0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))

df$species <- name
message(paste0("Finishing ",name, "!"))
#write.table with the df, include "species" as column
write.table(x = df, file = paste0(name,".mC_genes.tsv"), quote = F, sep = "\t")

######################################################################
## Read outputs from previous script and plot boxplot mCG/TPM, mCA/TPMs
######################################################################

list_of_files <- list.files(pattern = "mC_genes.tsv$",path = "~/Genes_to_CG_CA/", full.names = TRUE) %>% str_replace(pattern = "\\./", replacement = "")

all_species_gene_mC <- lapply(list_of_files, FUN = fread) %>% do.call(rbind,.)

all_species_gene_mC$decile[is.na(all_species_gene_mC$decile)] <- "No expr." 

all_species_gene_mC$decile <- factor(all_species_gene_mC$decile, levels = c("No expr.","0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))

all_species_gene_mC$species <- factor(all_species_gene_mC$species, levels = c("Honeybee","Octopus","Lancelet_neural","Lamprey","Shark","Zebrafish_forebrain","Xenopus","Gallus","Parus","Platypus","Opossum","Mouse","Human"))

cc <- scales::seq_gradient_pal("lightcyan", "royalblue4", "Lab")(seq(0,1,length.out=11))


gg_mCA_genes <- ggplot(all_species_gene_mC, aes(x = species, y = mCA, fill = decile)) + geom_boxplot(outlier.shape = NA) + 
  theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90),legend.position="none") + ylim(c(0,0.1)) + scale_fill_manual(values = cc)


gg_mCG_genes <- ggplot(all_species_gene_mC, aes(x = species, y = mCG, fill = decile)) + geom_boxplot(outlier.shape = NA) + 
  theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90),legend.position="none") + ylim(c(0,1)) + scale_fill_manual(values = cc)

gg_genes_combo <- plot_grid( gg_mCA_genes, gg_mCG_genes, gg_mCG_mock, ncol = 1 )

ggsave(filename = "BoxPlot_TPM_vs_mC_amphioxus.pdf", plot = gg_genes_combo, height = 10, width = 7)

# obtain spearman correlation coeficients for each species, in CG and CH context

for (i in unique(all_species_gene_mC$species)){
  all_species_gene_mC[all_species_gene_mC$species %in% c(i),] -> test_df
  message(i)
  cor(test_df$tpm, test_df$mCG, method = "spearman", use="pairwise") %>% message()
  cor(test_df$tpm, test_df$mCA, method = "spearman", use="pairwise") %>% message()
  cor(test_df$mCG, test_df$mCA, method = "spearman", use="pairwise") %>% message()
}


######################################################################
## Classify genes by mCA and mCG gene body level
######################################################################

compare_mC_deciles <- function(name){
  df <- fread(paste0("~/Genes_to_CG_CA/",name,".mC_genes.tsv"))
  df$mCA_decile <- cut(df$mCA, breaks = c(quantile(df$mCA, seq(0, 1, length = 11), type = 5)),
                       labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))
  df$mCA_decile[is.na(df$mCA_decile)] <- "0-10%" 
  df$mCG_decile <- cut(df$mCG, breaks = c(quantile(df$mCG, seq(0, 1, length = 11), type = 5)),
                       labels = c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))
  df$mCG_decile[is.na(df$mCG_decile)] <- "0-10%"
  write.table(x = df, file = paste0("~/Genes_to_CG_CA/",name,".mC_genes.deciles.tsv"), row.names = F, quote = F, sep = "\t")
  
  df <- df[!(is.na(df$mCG)) && !(is.na(df$mCA)) , ]
  
}

lapply(c("Human","Mouse","Opossum","Platypus","Parus","Gallus","Xenopus","Zebrafish_forebrain","Shark","Lamprey","Lancelet_neural","Octopus","Honeybee"), FUN = compare_mC_deciles)


######################################################################
## Obtain GO with gProfileR and plot enrichments
######################################################################

library(gProfileR)

gprofiler_species_specific <- function(name, species){
  dir <- "/Species_GOs/"
  df <- fread(paste0("~/Genes_to_CG_CA/",name,".mC_genes.deciles.tsv"))
  
  mCA_hi <- unique(df$target_id[df$mCA_decile == "90-100%"])
  mCA_low <- unique(df$target_id[df$mCA_decile == "0-10%"])
  mCG_hi <- unique(df$target_id[df$mCG_decile == "90-100%"])
  mCG_low <- unique(df$target_id[df$mCG_decile == "0-10%"])
  
  gprof_mCA_hi <- gprofiler(mCA_hi, organism = species, correction_method = "gSCS")
  gprof_mCA_low <- gprofiler(mCA_low, organism = species, correction_method = "gSCS")
  gprof_mCG_hi <- gprofiler(mCG_hi, organism = species, correction_method = "gSCS")
  gprof_mCG_low <- gprofiler(mCG_low, organism = species, correction_method = "gSCS")
  gprof_mCA_hi$type <- "mCA_hi" 
  gprof_mCA_low$type <- "mCA_low" 
  gprof_mCG_hi$type <- "mCG_hi" 
  gprof_mCG_low$type <- "mCG_low" 
  
  gos_all <- rbind(gprof_mCA_hi, gprof_mCA_low, gprof_mCG_hi, gprof_mCG_low)
  gos_all$species <- name
  
  write.table(x = gos_all, file = paste0(dir,name,"_gProfiler.species_specific.tsv"), quote = F, row.names = F, col.names = T, sep = "\t")
  #return(gos_all)
}

gprofiler_species_specific("Human", species = "hsapiens")
gprofiler_species_specific("Opossum", species = "mdomestica")
gprofiler_species_specific("Gallus", species = "ggallus")
gprofiler_species_specific("Mouse", species = "mmusculus")
gprofiler_species_specific("Platypus", species = "oanatinus")
gprofiler_species_specific("Xenopus", species = "xtropicalis")
gprofiler_species_specific("Octopus", species = "obimaculoides")
gprofiler_species_specific("Honeybee", species = "amellifera")
gprofiler_species_specific(name = "Parus", species = "pmajor")
gprofiler_species_specific(name = "Zebrafish_forebrain", species = "drerio")

import_GOs_ss_gprofiler <- function(name){
  file <- paste0("~/Species_GOs/",name,"_gProfiler.species_specific.tsv")
  df <- fread(file, sep = "\t") %>% .[,c("domain", "term.name","term.id","p.value","type","species")] %>%
    mutate(negative_log10_of_adjusted_p_value = -log10(p.value)) %>% dplyr::rename(term_name = term.name, term_id = term.id, adjusted_p_value = p.value,
                                                                                   source = domain)
  return(df)
}

all_GOs_ss <- lapply(c("Human","Mouse","Opossum","Platypus","Parus","Gallus","Xenopus","Zebrafish_forebrain","Octopus","Honeybee"), FUN = import_GOs_ss_gprofiler) %>% do.call(rbind,.) %>%
  data.frame()

import_GOs_orthofinder_gprofiler <- function(name){
  file <- paste0("~/Species_GOs/",name,"_gProfiler.orthofinder.tsv")
  df <- fread(file, sep = "\t") %>% .[,c("domain", "term.name","term.id","p.value","type","species")] %>%
    mutate(negative_log10_of_adjusted_p_value = -log10(p.value)) %>% dplyr::rename(term_name = term.name, term_id = term.id, adjusted_p_value = p.value,
                                                                                   source = domain)
  return(df)
}

all_GOs_of <- lapply(c("Shark","Lamprey"), FUN = import_GOs_orthofinder_gprofiler) %>% do.call(rbind,.) %>%
  data.frame()

all_GOs <- rbind(all_GOs_ss, all_GOs_of, amphi_neural_GOs)


#make table of all GOs and KEGGs

simplified_all_GOs <- all_GOs %>% dplyr::filter(!source %in% c("hp","hpa","mir","tf","HP","cor","rea")) %>% 
  dplyr::filter(adjusted_p_value < 0.05) %>%
  dplyr::mutate(source = str_replace(source,pattern = "GO:", replacement = "")) %>%
  .[,c(6,5,3,1,2,4,7)] %>% dplyr::rename(Gene_group = type)
write.table(x = simplified_all_GOs, row.names = F,
            file = "~/Dropbox/BrainZoo/Data/all_GOs_simplified.csv", quote = T, sep = ",")

# check most conserved GOs in the dataset for top mCA, top mCG, bottom mCA, bottom mCG (found at least in 4 species)
table(all_GOs$term_name[grepl(all_GOs$term_id, pattern = "GO") & grepl(all_GOs$type, pattern = "mCG_hi")]) %>% data.frame() %>% filter(Freq > 0) %>% arrange(Freq)
table(all_GOs$term_name[grepl(all_GOs$term_id, pattern = "GO") & grepl(all_GOs$type, pattern = "mCG_lo")]) %>% data.frame() %>% filter(Freq > 4) %>% arrange(Freq)
table(all_GOs$term_name[grepl(all_GOs$term_id, pattern = "GO") & grepl(all_GOs$type, pattern = "mCA_hi")]) %>% data.frame() %>% filter(Freq > 4) %>% arrange(Freq)
table(all_GOs$term_name[grepl(all_GOs$term_id, pattern = "GO") & grepl(all_GOs$type, pattern = "mCA_lo")]) %>% data.frame() %>% filter(Freq > 4) %>% arrange(Freq)

# import GOs to plot for the CA context and plot enrichments
gos_to_check <- fread(file = "conserved_mCA_GOs", sep = "\t", header = F) %>% .$V1

GOs_to_plot <- all_GOs[all_GOs$term_name %in% gos_to_check, ]

GOs_to_plot <- GOs_to_plot[!grepl(GOs_to_plot$term_id, pattern = "REAC"), ]

GOs_to_plot$long_name <- paste0(GOs_to_plot$term_name," (",GOs_to_plot$term_id, " - ", GOs_to_plot$source, ")" )
GOs_to_plot$species <- factor(GOs_to_plot$species, levels = c("Honeybee","Octopus","Lancelet_neural","Lamprey","Shark","Zebrafish_forebrain","Xenopus","Gallus","Parus","Platypus","Opossum","Mouse","Human"))
GOs_to_plot$term_name <- factor(GOs_to_plot$term_name, levels = gos_to_check)
GOs_to_plot$negative_log10_of_adjusted_p_value <- as.numeric(GOs_to_plot$negative_log10_of_adjusted_p_value)
GOs_to_plot$`-log10(q-val)` <- GOs_to_plot$negative_log10_of_adjusted_p_value
rename_GOs <- GOs_to_plot[,c("long_name","term_name")] %>% unique()

gg_GO_enrichments2 <- ggplot(GOs_to_plot[GOs_to_plot$type %in% c("mCA_hi","mCA_low"),], aes(species, term_name)) + geom_point(aes(size=`-log10(q-val)`, color = type), alpha=0.75) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90) ) + 
  scale_size_area(max_size = 8,breaks=c(3,5,10,50,100)) + 
  xlab("") + ylab("") + 
  scale_y_discrete(breaks=rename_GOs$term_name,
                   labels=rename_GOs$term_name) +
  facet_grid(.~type) +
  scale_fill_discrete(name = "New Legend Title")


ggsave(filename = "GO_enrichments_ss.pdf", plot = gg_GO_enrichments2, height = 10, width = 10)

# import GOs to plot for the CG context and plot enrichments

all_GOs_CG <- all_GOs[!all_GOs$type %in% c("mCA_hi","mCA_low"),]

gos_to_check_CG <- fread(file = "mCG_GOs", sep = "\t", header = F) %>% .$V1

GOs_to_plot_CG <- all_GOs_CG[all_GOs_CG$term_name %in% gos_to_check_CG, ]

GOs_to_plot_CG <- GOs_to_plot_CG[!grepl(GOs_to_plot_CG$term_id, pattern = "REAC"), ]

GOs_to_plot_CG$long_name <- paste0(GOs_to_plot_CG$term_name," (",GOs_to_plot_CG$term_id, " - ", GOs_to_plot_CG$source, ")" )
GOs_to_plot_CG$species <- factor(GOs_to_plot_CG$species, levels = c("Honeybee","Octopus","Lancelet_neural","Lamprey","Shark","Zebrafish_forebrain","Xenopus","Gallus","Parus","Platypus","Opossum","Mouse","Human"))
GOs_to_plot_CG$term_name <- factor(GOs_to_plot_CG$term_name, levels = gos_to_check_CG)
GOs_to_plot_CG$negative_log10_of_adjusted_p_value <- as.numeric(GOs_to_plot_CG$negative_log10_of_adjusted_p_value)
GOs_to_plot_CG$`-log10(q-val)` <- GOs_to_plot_CG$negative_log10_of_adjusted_p_value
rename_GOs_CG <- GOs_to_plot_CG[,c("long_name","term_name")] %>% unique()

gg_GO_enrichments_CG <- ggplot(GOs_to_plot_CG, aes(species, term_name)) + geom_point(aes(size=`-log10(q-val)`, color = type), alpha=0.75) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90) ) + 
  scale_size_area(max_size = 8,breaks=c(3,5,10,50,100)) + 
  xlab("") + ylab("") + 
  scale_y_discrete(breaks=rename_GOs_CG$term_name,
                   labels=rename_GOs_CG$term_name) +
  facet_grid(.~type) +
  scale_fill_discrete(name = "New Legend Title")

ggsave(filename = "~/Dropbox/BrainZoo/Plots/GO_enrichments_CG_ss.pdf", plot = gg_GO_enrichments_CG, height = 10, width = 10)


######################################################################
## Overlaps per categories (top mCA vs bottom mCG, etc.)
######################################################################

list_of_files <- list.files(pattern = "mC_genes.deciles.tsv$",path = "~/Genes_to_CG_CA/", full.names = TRUE) %>% str_replace(pattern = "\\./", replacement = "")

all_species_gene_mC_deciles <- lapply(list_of_files, FUN = fread) %>% do.call(rbind,.)
all_species_gene_mC_deciles$species <- factor(all_species_gene_mC_deciles$species, levels = c("Honeybee","Octopus","Lancelet_neural","Lamprey","Shark","Zebrafish_forebrain","Xenopus","Gallus","Parus","Platypus","Opossum","Mouse","Human"))
species_order <- c("Honeybee","Octopus","Lancelet_neural","Lamprey","Shark","Zebrafish_forebrain","Xenopus","Gallus","Parus","Platypus","Opossum","Mouse","Human")

intersect_genes_by_deciles <- function(species, setCG = "0-10%", setCA = "90-100%"){
  genesA <- all_species_gene_mC_deciles$target_id[all_species_gene_mC_deciles$mCG_decile == setCG & all_species_gene_mC_deciles$species == species]
  genesB <- all_species_gene_mC_deciles$target_id[all_species_gene_mC_deciles$mCA_decile == setCA & all_species_gene_mC_deciles$species == species]
  both <- sum(genesA %in% genesB)
  A_uni <- sum(!genesA %in% genesB)
  B_uni <- sum(!genesB %in% genesA)
  df <- data.frame(value = c(both, A_uni, B_uni), type = c("overlap",paste0("mCG_",setCG),paste0("mCA_",setCA) )) %>%
    mutate(perc = 100*value/sum(value))
  df$type <- factor(df$type, levels = c(paste0("mCG_",setCG),"overlap",paste0("mCA_",setCA)))
  df$species <- species
  return(df)
}

hiCA_lowCG <- lapply( species_order, intersect_genes_by_deciles) %>% do.call(rbind,.) %>% data.frame()
hiCA_lowCG$species <- factor(hiCA_lowCG$species, levels = species_order)
hiCA_hiCG <- lapply( species_order, intersect_genes_by_deciles, setCG = "90-100%", setCA = "90-100%") %>% do.call(rbind,.) %>% data.frame()
hiCA_hiCG$species <- factor(hiCA_hiCG$species, levels = species_order)
lowCA_lowCG <- lapply( species_order, intersect_genes_by_deciles, setCG = "0-10%", setCA = "0-10%") %>% do.call(rbind,.) %>% data.frame()
lowCA_lowCG$species <- factor(lowCA_lowCG$species, levels = species_order)
lowCA_highCG <- lapply( species_order, intersect_genes_by_deciles, setCG = "90-100%", setCA = "0-10%") %>% do.call(rbind,.) %>% data.frame()
lowCA_highCG$species <- factor(lowCA_highCG$species, levels = species_order)

cbPalette <- c( "#0072B2", "#009E73", "#D55E00")


plot_overlaps <- function(x){
  ggplot(x[x$type == "overlap",], aes(x = species, y = perc)) + geom_bar( stat = "identity" ) + 
    scale_fill_manual(values = cbPalette) +
    theme( axis.text.x = element_text(angle = 90, hjust = 0)
           ,plot.background = element_blank()
           ,panel.grid.major = element_blank()
           ,panel.grid.minor = element_blank()
           ,panel.border = element_blank()
    ) +
    # X axis label
    xlab(label = "") + 
    # Y axis label
    ylab(label = "percentage (%)") +
    ylim(c(0,35))
}

overlap_deciles_plot <- plot_grid( plot_overlaps(hiCA_lowCG),
                                   plot_overlaps(hiCA_hiCG),
                                   plot_overlaps(lowCA_lowCG),
                                   plot_overlaps(lowCA_highCG), ncol = 2)

ggsave(overlap_deciles_plot, filename = "Overlap_mC_deciles.pdf")


