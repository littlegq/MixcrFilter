library(ggplot2)
library(scales)
library(base)
library(grid)
library(gridExtra)
library(Rmisc)
library(reshape2)

# List of Input files:
# rpm_ex17_51bp_AITLs.txt
# combine.AB.top10.CDR3reads100.flt.txt
# SampleID.index.txt
# The whole flt.clones folder
# TCR138.GRCh38.fpkm.txt
# rpm.txt
# clonal_check.txt

# Exclude error cases
ex_cases <- c("NKCL.Chan-12", "NKCL.Chan-53", "pNOS-NOS1", "ALCL-Cellline.JB6",
              "AITL.Yoo.PAT13.RNA", "pNOS-S93-6748", "pNOS-S98-1772")
## Read dictionary files for case IDs used for publication
pubid <- read.table("SampleID.index.txt", header = T, stringsAsFactors = F)
pubid <- pubid[which(!pubid$Conversed %in% ex_cases), ]

# Determine the case IDs with low TCR expression
## Read data from the file
d17 <- read.table("rpm_ex17_51bp_AITLs.txt", header = T)

## exclude known cases with errors
d17 <- d17[which(!d17$Sample %in% ex_cases, d17$Type != "09.gdTL"),]  
rpm <- melt(d17[, c('TRA.rpm', 'TRB.rpm', 'TRG.rpm', 'TRD.rpm', 'Type')], id = 'Type')
colnames(rpm)[2:3] <- c('strand', 'RPM')

## Get the cutoff from the NKCL cases 
nkcl_rpm <- rpm[which(rpm$Type == "07.NKCL"),]
nkcl_TRA_maxRPM <- max(nkcl_rpm[which(nkcl_rpm$strand == "TRA.rpm"), 'RPM'])
nkcl_TRB_maxRPM <- max(nkcl_rpm[which(nkcl_rpm$strand == "TRB.rpm"), 'RPM'])

## IDs of excluded cases
lowTRA <- d17$TRA.rpm <= nkcl_TRA_maxRPM
lowTRB <- d17$TRB.rpm <= nkcl_TRB_maxRPM
nkcl_TRA_maxRPM
nkcl_TRB_maxRPM
lowTCRsampleIDs <- d17[lowTRA & lowTRB, 'SampleID']
TRA.excluded <- d17[lowTRA & !lowTRB, 'SampleID']
TRB.excluded <- d17[!lowTRA & lowTRB, 'SampleID']

# Figure 2: Cutoff figure
## Read the ID table
d <- merge(d17, pubid, by.x = 'SampleID', by.y = 'Conversed')

## Extract necessary information of the barplot
d <- d[, c('TRA.rpm', 'TRB.rpm', 'Type', 'Publish')]
colnames(d)[4] <- 'SampleID'

## Plot TRA
idorder <- d[order(d$Type, d$TRA.rpm), ]$SampleID
pa <- ggplot(d, aes(x = SampleID, y = TRA.rpm, fill = Type))
pa <- pa + geom_bar(stat = "identity") +  
  scale_x_discrete(limits = idorder) +
  geom_hline(yintercept = nkcl_TRA_maxRPM) +
  theme_bw() + 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1)),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = rel(0.4))) 

## Plot TRB
idorder <- d[order(d$Type, d$TRB.rpm), ]$SampleID
pb <- ggplot(d,aes(x = SampleID, y = TRB.rpm, fill = Type))
pb <- pb + geom_bar(stat = "identity") +  
  scale_x_discrete(limits = idorder) +
  geom_hline(yintercept = nkcl_TRB_maxRPM) +
  theme_bw() + 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1)),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = rel(0.4)))

# lowTCRsampleIDs <- pubid[pubid$Conversed %in% lowTCRsampleIDs, 'Publish']

## Plot Figure 2
pdf("Figure2.cutoff.pdf", height = 8.5, width = 11)
p <- list(pa, pb)
lo <- matrix(
  c(1, 2), 
  nrow = 2, byrow = T)
multiplot(plotlist = p, layout = lo)
dev.off()

## clean up variables
rm(lo, nkcl_rpm, idorder, lowTRA, lowTRB, p, pa, pb)

#########################################################
# Classify mono-clonal and poly-clonal cases
min_ratio <- 10
min_top2 <- 0.1

cv <- read.table("TCR138.GRCh38.fpkm.txt", header = T)
TRA.FPKM1 <- c()
TRA.FPKM2 <- c()
TRA.FPKM3 <- c()
TRA.top2 <- c()
TRA.FPKM.p1 <- c()
TRA.FPKM.p2 <- c()
TRA.FPKM.p3 <- c()
TRB.FPKM1 <- c()
TRB.FPKM2 <- c()
TRB.FPKM3 <- c()
TRB.top2 <- c()
TRB.FPKM.p1 <- c()
TRB.FPKM.p2 <- c()
TRB.FPKM.p3 <- c()

for (id in pubid$TCRV) {
  ordv <- cv[which(cv$Type %in% c('Type', 'AV')),
             c('Gene', id)]
  ordv <- ordv[order(-ordv[, 2]),]
  TRA.FPKM1 <- c(TRA.FPKM1, ordv[1, 2])
  TRA.FPKM2 <- c(TRA.FPKM2, ordv[2, 2])
  TRA.FPKM3 <- c(TRA.FPKM3, ordv[3, 2])
  TRA.top2 <- c(TRA.top2, (ordv[1, 2] + ordv[2, 2]) / sum(ordv[, 2]))
  TRA.FPKM.p1 <- c(TRA.FPKM.p1, ordv[1, 2] / sum(ordv[, 2]))
  TRA.FPKM.p2 <- c(TRA.FPKM.p2, ordv[2, 2] / sum(ordv[, 2]))
  TRA.FPKM.p3 <- c(TRA.FPKM.p3, ordv[3, 2] / sum(ordv[, 2]))
  ordv <- cv[which(cv$Type %in% c('Type', 'BV')),
             c('Gene', id)]
  ordv <- ordv[order(-ordv[, 2]),]
  TRB.FPKM1 <- c(TRB.FPKM1, ordv[1, 2])
  TRB.FPKM2 <- c(TRB.FPKM2, ordv[2, 2])
  TRB.FPKM3 <- c(TRB.FPKM3, ordv[3, 2])
  TRB.top2 <- c(TRB.top2, (ordv[1, 2] + ordv[2, 2]) / sum(ordv[, 2]))
  TRB.FPKM.p1 <- c(TRB.FPKM.p1, ordv[1, 2] / sum(ordv[, 2]))
  TRB.FPKM.p2 <- c(TRB.FPKM.p2, ordv[2, 2] / sum(ordv[, 2]))
  TRB.FPKM.p3 <- c(TRB.FPKM.p3, ordv[3, 2] / sum(ordv[, 2]))
}
v_eval <- data.frame(
  pubid$Conversed, pubid$Publish,
  TRA.FPKM1, TRA.FPKM2, TRA.FPKM3, TRA.FPKM.p1, TRA.FPKM.p2, TRA.FPKM.p3, TRA.top2, 
  TRB.FPKM1, TRB.FPKM2, TRB.FPKM3, TRB.FPKM.p1, TRB.FPKM.p2, TRB.FPKM.p3, TRB.top2
)
v_eval$TRA.FPKM.stat <- 
  (v_eval$TRA.FPKM1 + v_eval$TRA.FPKM2 >= v_eval$TRA.FPKM3 * min_ratio) &
  (!v_eval$pubid.Conversed %in% TRA.excluded)
v_eval$TRB.FPKM.stat <- 
  (v_eval$TRB.FPKM1 + v_eval$TRB.FPKM2 >= v_eval$TRB.FPKM3 * min_ratio) &
  (!v_eval$pubid.Conversed %in% TRB.excluded)

## Exclude error cases
v_eval <- v_eval[!v_eval$pubid.Conversed %in% ex_cases, ]
v_eval <- v_eval[, c(2:18)]

## remove used variables
rm(cv, ordv, id,
   TRA.top2, TRA.FPKM1, TRA.FPKM2, TRA.FPKM3, TRA.FPKM.p1, TRA.FPKM.p2, TRA.FPKM.p3,
   TRB.top2, TRB.FPKM1, TRB.FPKM2, TRB.FPKM3, TRB.FPKM.p1, TRB.FPKM.p2, TRB.FPKM.p3)

## read clonality from CDR3 information
cdr3 <- read.table('clonal_check.txt')
colnames(cdr3) <- c(
  'Sample',
  'Strand1', 'TRA.stat', 'TRA.CDR3.1', 'TRA.CDR3.2', 'TRA.CDR3.3', 'TRA.p1', 'TRA.p2', 'TRA.p3',
  'Strand2', 'TRB.stat', 'TRB.CDR3.1', 'TRB.CDR3.2', 'TRB.CDR3.3', 'TRB.p1', 'TRB.p2', 'TRB.p3'
)
cdr3$TRA.CDR3.stat <- 
  ((cdr3$TRA.CDR3.1 + cdr3$TRA.CDR3.2 >= cdr3$TRA.CDR3.3 * min_ratio) | 
  (cdr3$TRA.p1 + cdr3$TRA.p2 >= min_top2)) &
  (!cdr3$Sample %in% TRA.excluded)
cdr3$TRB.CDR3.stat <- 
  ((cdr3$TRB.CDR3.1 + cdr3$TRB.CDR3.2 >= cdr3$TRB.CDR3.3 * min_ratio) |
  (cdr3$TRB.p1 + cdr3$TRB.p2 >= min_top2)) &
  (!cdr3$Sample %in% TRB.excluded)

## Exclude error cases
cdr3 <- cdr3[!cdr3$Sample %in% ex_cases,]
cdr3 <- merge(cdr3,pubid,by.x = 'Sample', by.y = 'Conversed')
cdr3 <- cdr3[, c(2, 4:10, 12:19, 21)]

# Compare the two evaluation
## Determine the pubID of low TCR cases
lowTCR_pubid <- pubid[pubid$Conversed %in% lowTCRsampleIDs, 'Publish']

## Merge the evaluation data with CDR3
stat_compare <- merge(v_eval, cdr3, by.x = 'pubid.Publish', by.y = 'Publish')

## Exclude cases with low TCR expression
stat_compare <- stat_compare[!stat_compare$pubid.Publish %in% lowTCR_pubid,]

## Make the comparison between CDR3 and FPKM systems
attach(stat_compare)
stat_compare$TRA.con <- TRA.FPKM.stat == TRA.CDR3.stat
stat_compare$TRB.con <- TRB.FPKM.stat == TRB.CDR3.stat
stat_compare$cons <- (TRA.FPKM.stat | TRB.FPKM.stat) == (TRA.CDR3.stat | TRB.CDR3.stat)
detach(stat_compare)

## Output the comparing result
write.table(stat_compare, 'tableS2.TCR.clonal.compare.txt', row.names = F, sep = "\t")

## Table with identifiable IDs
stat_compare_oriID <- merge(stat_compare, pubid, by.x = "pubid.Publish", 
                            by.y = "Publish")
stat_compare_oriID <- stat_compare_oriID[, c(37, 2:36)]
write.csv(stat_compare_oriID, 'tableS1.TCR.clonal.compare.oriID.csv', row.names = F)

# Determine case orders of the plots
stat_compare$monoclonal <- with(stat_compare, TRA.CDR3.stat | TRB.CDR3.stat)
stat_compare$CDR3.mean.top2 <- 
  with(stat_compare, (TRA.p1 + TRA.p2 + TRB.p1 + TRB.p2) / 2)
TRA.excluded <- pubid[pubid$Conversed %in% TRA.excluded, 'Publish']
TRB.excluded <- pubid[pubid$Conversed %in% TRB.excluded, 'Publish']

## calculate mean CDR3 top2 clone size regardless of TCR negative strands
neg_select <- stat_compare$pubid.Publish %in% TRA.excluded
stat_compare[neg_select, 'CDR3.mean.top2'] <- with(stat_compare[neg_select,], TRB.p1 + TRB.p2)
neg_select <- stat_compare$pubid.Publish %in% TRB.excluded
stat_compare[neg_select,'CDR3.mean.top2'] <- with(stat_compare[neg_select,], TRA.p1 + TRA.p2)
rm(neg_select)

## Get the clonal proportion of AITL from the CDR3 results
aitl_stat <- stat_compare[grep("AITL", stat_compare$pubid.Publish), ]
aitl_stat <- summary(aitl_stat$monoclonal)

summary(stat_compare$TRA.con)
# Mode   FALSE    TRUE    NA's 
# logical      23      64       0
summary(stat_compare$TRB.con)
# Mode   FALSE    TRUE    NA's 
# logical      13      74       0 

## Determine the case orders of the plots
stat_compare[stat_compare$pubid.Publish %in% 
               paste0("T-cell-", 1:4), 'monoclonal'] <- "Control-1"
stat_compare[stat_compare$pubid.Publish %in% 
               paste0("T-cell-", 5:7), 'monoclonal'] <- "Control-2"
stat_compare[stat_compare$pubid.Publish %in% 
               paste0("T-cell-", 8:10), 'monoclonal'] <- "Control-3"
case_orders <- stat_compare[
  order(stat_compare$monoclonal, stat_compare$CDR3.mean.top2), 'pubid.Publish']

## Clean up used variables
rm(v_eval, cdr3, rpm)

###############################################
# Figure 3
## Read data from the file
d <- read.table("combine.AB.top10.CDR3reads100.flt.txt")
colnames(d) <- c('OldID', 'Sample', 'Strand', 'Order', 'Count', 'Portion')

## Define case types
d$Type[grep("AITL", d$Sample)] <- 'T1.AITL'
d$Type[grep("ALCL", d$Sample)] <- 'T2.ALCL'
d$Type[grep("pNOS", d$Sample)] <- 'T4.pNOS'
d$Type[grep("NKCL", d$Sample)] <- 'NK1.NKCL'
d$Type[grep("NormalNK", d$Sample)] <- 'NK0.NormalNK'
d$Type[grep("NK.Cellline", d$Sample)] <- 'NK2.Cellline'
d$Type[grep("ALCL-Cellline", d$Sample)] <- 'T3.ALCL-Cellline'
d$Type[grep("gdTL", d$Sample)] <- 'T5.gdTL'
d$Type[grep("NAIVE", d$Sample)] <- 'T01.NAIVE'
d$Type[grep("TFH", d$Sample)] <- 'T02.TFH'
d$Type[grep("TEFF", d$Sample)] <- 'T03.TEFF'
d$Type <- as.factor(d$Type)
d$TotalCount <- round(d$Count / d$Portion)

## Exclude cases with errors
d <- d[!d$Sample %in% ex_cases,]

## Exclude the cases with insufficient TCR expression
d <- d[(!d$Sample %in% lowTCRsampleIDs), ]

## Define the color scale
colors = c("1" = "darkred",
           "2" = "red",
           "3" = "orange",
           "4" = "yellow",
           "5" = "lightgreen",
           "6" = "green",
           "7" = "lightblue",
           "8" = "blue",
           "9" = "darkblue",
           "10" = "purple",
           "11" = "grey",
           "0" = "white")

## Prepare data for the TRB plot
i <- grep("AITL|0_", d$Sample)
b <- subset(d[i, ], Strand == "TRB")
b <- merge(b, pubid, by.x = 'Sample', by.y = 'Conversed')
colnames(b)[1] <- 'Conversed_ID'
colnames(b)[10] <- 'Sample'
b[b$Sample %in% TRB.excluded, 'Order'] <- 0   # fill the low express ones with white
b <- within(b, Sample <- factor(Sample, levels = case_orders)) 

## Plot TRB
pb <- ggplot(b, aes(x = factor(Sample), y = Portion, fill = factor(Order), width = .5))
pb <- pb + geom_bar(stat = "identity", color = "black") + scale_y_continuous(labels = percent)
pb <- pb + scale_fill_manual(values = colors)
pb <- pb + theme_bw() + theme(
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
  legend.position = "none",
  panel.grid.major.y = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.grid.minor.x = element_blank()) 

## Prepare data for the TRA plot
a <- subset(d[i, ], Strand == "TRA")
a[a$Sample %in% TRA.excluded, 'Order'] <- 0   # fill the low express ones with white
a <- merge(a, pubid, by.x = 'Sample', by.y = 'Conversed')
colnames(a)[1] <- 'Conversed_ID'
colnames(a)[10] <- 'Sample'
a[a$Sample %in% TRA.excluded, 'Order'] <- 0   # fill the low express ones with white
a <- within(a, Sample <- factor(Sample, levels = case_orders))  

## Plot TRA
pa <- ggplot(a, aes(x = factor(Sample), y = Portion, fill = factor(Order), width = .5))
pa <- pa + geom_bar(stat = "identity", color = "black") + scale_y_continuous(labels = percent)
pa <- pa + scale_fill_manual(values = colors)
pa <- pa + theme_bw() + theme(
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
  legend.position = "none",
  panel.grid.major.y = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.grid.minor.x = element_blank()) 

## Case IDs for the lower panels (CDR3)
ids <- c('AITL-6-592-4180',
         'AITL-Ly19',
         'Yoo.PAT6.RNA')

## Read the data from the flt.clones folder and plot the CDR3 barplots
class <- 'flt'
p <- list()
pidx <- 0
for (gene in c("TRA", "TRB")){
  for (id in ids){
    cdr3 <- read.table(
      paste(class, ".clones/", gene, "/", id, ".", gene, ".", class, ".clones.txt", sep = ""),
      header = T)
    if(nrow(cdr3) == 0) next
    cdr3 <- cdr3[1:10, ]
    colnames(cdr3)[8] <- 'CDR3.AA'
    cdr3 <- within(cdr3,CDR3.AA <- factor(CDR3.AA, levels = CDR3.AA))
    pidx <- pidx + 1
    p[[pidx]] <- ggplot(
      cdr3, aes(x = as.factor(as.numeric(rownames(cdr3))),
               y = Clone_count)) +
      geom_bar(stat = "identity", width = .7) + 
      theme_bw() +
      theme(panel.grid.major.y = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())   
  }
}

# Read FPKM data
cv <- read.table("TCR138.GRCh38.fpkm.txt", header = T)
for (gene in c("AV", "BV")){
  for (id in pubid[which(pubid$Original %in% ids), 'TCRV']){
    ordv <- cv[which(cv$Type %in% c('Type', gene)), c('Gene', id)]
    ordv <- ordv[order(-ordv[, 2]), ]
    ordv <- ordv[1:10, ]
    colnames(ordv) <- c('TCRV', 'FPKM')
    rownames(ordv) <- c(1:10)
    pidx <- pidx + 1
    p[[pidx]] <- ggplot(
      ordv, aes(x = as.factor(as.numeric(rownames(ordv))), y = FPKM)) +
      geom_bar(stat = "identity",width=.7) + 
      theme_bw() +
      theme(panel.grid.major.y = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())   
  }
}

df <- data.frame(
  variable = c("Monal-clonal", "Poly-clonal"),
  value = as.numeric(aitl_stat[2:3])
)

## Pie chart
pc <- ggplot(df, aes(x = "", y = value, fill = variable)) +
  geom_bar(width = 1, stat = "identity", col = "black") +
  scale_fill_manual(values = c("grey", "white")) +
  coord_polar("y", start = pi / 3) + 
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") 

## Plot Figure 3
p[[13]] <- pa
p[[14]] <- pb
p[[15]] <- pc
pdf("Figure3.pdf", height = 22, width = 17)
lo <- matrix(
  c(13, 13, 13, 13, 13,
    13, 13, 13, 13, 13,
    14, 14, 14, 14, 14,
    14, 14, 14, 14, 14,
    1, 4, 7, 10, 15, 
    2, 5, 8, 11, 15,
    3, 6, 9, 12, 15), 
  nrow = 7, byrow = T)
multiplot(plotlist = p, layout = lo)
dev.off()

## Clean up used variables
rm(a, b, cv, df, lo, ordv, aitl_stat, class, gene, i, id, ids, p, pa, pb, pc, pidx)


###############################################
# Figure 6 (ALCL)
## Extract ALCL index from the dataframe
i <- grep("ALCL|0_",d$Sample)

## Plot TRB
b <- subset(d[i,], Strand == "TRB")
b <- merge(b, pubid, by.x = 'Sample', by.y = 'Conversed')
colnames(b)[1] <- 'Conversed_ID'
colnames(b)[10] <- 'Sample'
b[b$Sample %in% TRB.excluded, 'Order'] <- 0   # fill the low express ones with white
b <- within(b, Sample <- factor(Sample, levels = case_orders))  
pb <- ggplot(b, aes(x = factor(Sample), y = Portion, fill = factor(Order), width = .5))
pb <- pb + geom_bar(stat = "identity", color = "black") + scale_y_continuous(labels = percent)
pb <- pb + scale_fill_manual(values = colors)
pb <- pb + theme_bw() + theme(
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
  legend.position = "none",
  panel.grid.major.y = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.grid.minor.x = element_blank()) 

## Plot TRA
a <- subset(d[i, ], Strand == "TRA")
a <- merge(a, pubid, by.x = 'Sample', by.y = 'Conversed')
colnames(a)[1] <- 'Conversed_ID'
colnames(a)[10] <- 'Sample'
a[a$Sample %in% TRA.excluded, 'Order'] <- 0   # fill the low express ones with white
a <- within(a, Sample <- factor(Sample, levels = case_orders))  
pa <- ggplot(a, aes(x = factor(Sample), y = Portion, fill = factor(Order), width = .5))
pa <- pa + geom_bar(stat = "identity", color = "black") + scale_y_continuous(labels = percent)
pa <- pa + scale_fill_manual(values = colors)
pa <- pa + theme_bw() + theme(
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
  legend.position = "none",
  panel.grid.major.y = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.grid.minor.x = element_blank()) 

## Prepare data for FPKM boxplot
df <- read.table("rpm.txt", header = T)

## Exclude cases with errors
df <- df[!df$SampleID %in% ex_cases,]
fpkm <- melt(df[, c('CD247.fpkm', 'CD3D.fpkm', 'CD3E.fpkm', 'CD3G.fpkm', 'Type')], id = 'Type')
colnames(fpkm)[2:3] <- c('gene', 'FPKM')

## Prepare data for max FPKM of V genes boxplot
vmax <- melt(d17[, c('TRAV.fpkm', 'TRBV.fpkm', 'TRGV.fpkm', 'TRDV.fpkm', 'Type')], id = 'Type')
colnames(vmax)[2:3] <- c('strand', 'FPKM')

## Define the subtype names for the plot
subtypes <- c('01.Normal_T' = 'Normal T cells',
              '02.AITL'   = 'AITL',
              '03.PTCL-NOS' = 'PTCL-NOS',
              '04.ALCL' = 'ALCL',
              '04.ALCL.noTCR' = 'ALCL(no TCR)',
              '04.ALCL.TCR' = 'ALCL(TCR)',
              '05.ALCL_Cellline' = 'ALCL Cell lines',
              '06.Normal_NK' = 'Normal NK cells',
              '07.NKCL' = 'NKCL',
              '08.NK_Cell_line' = 'NKCL Cell lines',
              '09.gdTL' = bquote(paste(gamma,delta,'-TCL')))
cd3genes <- list(bquote(paste('CD3', zeta)), 
                 bquote(paste('CD3', delta)), 
                 bquote(paste('CD3', epsilon)),
                 bquote(paste('CD3', gamma)))
strands <- list(bquote(paste('TCR', alpha)), 
                bquote(paste('TCR', beta)), 
                bquote(paste('TCR', gamma)),
                bquote(paste('TCR', delta)))

## Generate the boxplot for FPKM
pr <- ggplot(vmax, aes(x = Type, y = FPKM, color = strand))
pr <- pr + geom_boxplot(outlier.size = 1) + 
  ylim(0, quantile(vmax$FPKM, probs = 0.99)) + 
  theme_bw() + 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1)),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(labels = subtypes) +
  scale_colour_manual(values = 2:5, labels = strands)

## Generate the boxplot for FPKM
pf <- ggplot(fpkm, aes(x = Type, y = FPKM, color = gene))
pf <- pf + geom_boxplot(outlier.size = 1) + 
  ylim(0, quantile(fpkm$FPKM, probs = 0.995)) + 
  theme_bw() + 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1)),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(labels = subtypes) +
  scale_colour_manual(values = 2:5, labels = cd3genes)

pdf("Figure6.pdf", height = 22, width = 17)
p <- list(pa, pb, pr, pf)
lo <- matrix(
  c(1, 1,
    2, 2,
    3, 4), 
  nrow = 3, byrow = T)
multiplot(plotlist = p, layout = lo)
dev.off()

## Clean up used variables
rm(a, b, cdr3, df, fpkm_exLowTCR, lo, ALCL.lowTCR, ALCL.TCR, cd3genes, i, 
   p, pa, pb, pf, pr, strands, subtypes)

########################################################
# Figure 4 (PTCL-NOS)
## Extract case index from the dataframe
i <- grep("pNOS|0_", d$Sample)

## Prepare data for TRB barplot
b <- subset(d[i, ], Strand == "TRB")
b <- merge(b, pubid, by.x = 'Sample', by.y = 'Conversed')
colnames(b)[1] <- 'Conversed_ID'
colnames(b)[10] <- 'Sample'
b[b$Sample %in% TRB.excluded, 'Order'] <- 0   # fill the low express ones with white
b <- within(b, Sample <- factor(Sample, levels = case_orders))  

## Plot TRB
pb <- ggplot(b, aes(x = factor(Sample), y = Portion, fill = factor(Order), width=.5))
pb <- pb + geom_bar(stat = "identity", color = "black") + scale_y_continuous(labels = percent)
pb <- pb + scale_fill_manual(values = colors)
pb <- pb + theme_bw() + theme(
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
  legend.position = "none",
  panel.grid.major.y = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.grid.minor.x = element_blank()) 

## Prepare data from TRA barplot
a <- subset(d[i, ], Strand == "TRA")
a <- merge(a, pubid, by.x = 'Sample', by.y = 'Conversed')
colnames(a)[1] <- 'Conversed_ID'
colnames(a)[10] <- 'Sample'
a[a$Sample %in% TRA.excluded, 'Order'] <- 0   # fill the low express ones with white
a <- within(a, Sample <- factor(Sample, levels = case_orders))

## Plot TRA
pa <- ggplot(a, aes(x = factor(Sample), y = Portion, fill = factor(Order), width=.5))
pa <- pa + geom_bar(stat = "identity", color = "black") + scale_y_continuous(labels = percent)
pa <- pa + scale_fill_manual(values = colors)
pa <- pa + theme_bw() + theme(
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
  legend.position = "none",
  panel.grid.major.y = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.grid.minor.x = element_blank()) 

## Plot Figure 4
pdf("Figure4.PTCL-NOS.pdf",width = 11,height = 8.5)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(pa, vp = vplayout(1, 1))  # key is to define vplayout
print(pb, vp = vplayout(2, 1))
dev.off()

## Clean up used variables
rm(a, b, i, pa, pb, vplayout)

## Plot Figure 1 (Bi-allelic cases)
dt3 <- read.table("TableS3.txt", header = T)
library(ggplot2)
c12yes <- dt3[which(dt3$OOF == "Yes"), 'C12ratio']
c12no <- dt3[which(dt3$OOF == "No"), 'C12ratio']
t.test(c12yes, c12no, alternative = "greater")$p.value # 0.0008202897
p <- ggplot(dt3, aes(OOF, C12ratio))
p <- p + geom_boxplot() + geom_jitter(shape = 19) + theme_bw() + theme(
  panel.grid.major.y = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.grid.minor.x = element_blank()) 
ggsave(plot = p, height = 6, width = 6, dpi = 300, filename = "Figure1.NMD.pdf", useDingbats = FALSE)

## Test significance of the mutation frequency
mutFreq <-
  matrix(c(31, 5, 1, 3),
         nrow = 2,
         dimnames = list(Monoclonal = c("Mutated", "WildType"),
                         Polyclonal = c("Mutated", "WildType")))
fisher.test(mutFreq)$p.value # 0.02037422
fisher.test(mutFreq, alternative = "greater")$p.value # 0.02037422

## justify the ratio cutoff
d <- read.table("tableS2.txt", header = T)
d$TotalA <- d$TRA.CDR3.1 / d$TRA.p1
d$TotalB <- d$TRB.CDR3.1 / d$TRB.p1
total <- c(d$TotalA, d$TotalB)
d$TRAratio <- (d$TRA.CDR3.1 + d$TRA.CDR3.2) / d$TRA.CDR3.3
d$TRBratio <- (d$TRB.CDR3.1 + d$TRB.CDR3.2) / d$TRB.CDR3.3
for (i in 1:nrow(d)){
  d$maxRatio[i] <- max(d$TRAratio[i], d$TRBratio[i])
} 
d$Clonality <- (d$TRA.p1 + d$TRA.p2 >= 0.1) | (d$TRB.p1 + d$TRB.p2 >= 0.1)
ratio <- d[d$Clonality, 'maxRatio']
d$Clonality <- d$TRAratio >= 10 | d$TRBratio >= 10 | 
  (d$TRA.p1 + d$TRA.p2 >= 0.1) | (d$TRB.p1 + d$TRB.p2 >= 0.1)
x <- d[d$Clonality & d$TRAratio < 10 & d$TRBratio < 10, ]

d[which(d$Type == "normalT"), 'Clonality'] <- 'Normal-T'
p1 <- ggplot(d, aes(TRAratio, TRBratio, color = Clonality))
p1 <- p1 + geom_point() + scale_x_continuous(trans = "log10") + scale_y_continuous(trans = "log10")
for (i in c(6, 10, 14)) {
  h <- data.frame(x1 = 0, x2 = i, y1 = i, y2 = i)
  v <- data.frame(x1 = i, x2 = i, y1 = 0, y2 = i)
  p1 <- p1 + 
    geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = h) +
    geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = v)
}
ggsave(plot = p1, height = 8.5 * 0.8, width = 11 * 0.8, dpi = 300, 
       filename = "RatioCutoff.pdf", useDingbats = FALSE)
