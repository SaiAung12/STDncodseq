### excluding both Chloroplast and Mitochondria
library(phyloseq)
library(qiime2R)

# commands to load libraries and datasets
library(knitr)
library(ALDEx2)
library(CoDaSeq)
library(zCompositions)
library(igraph)
library(car)
library(grDevices)
library(propr)
library(vegan)
library(tidyr)
library(dplyr)
library(compositions)
library(xlsx)

##### STDnCoda_pipeline

#### addition for standard
library(ggplot2)
library(DESeq2)
library(dendextend)
library(microbiome) # data analysis and visualisation
library(microbiomeutilities) # some utility tools
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(DT) # interactive tables in html and markdown
library(data.table) # alternative to data.frame
library(ggsci)
library(scales)
library(DT)
library(htmlwidgets)

### READING IN OUR DATA

### from from STDncodseq03-01
D4STD0301 <- new.env(parent = emptyenv())

load("/home/pvannameihome/PG0xxx03-01/dada2fxxx/xxx_aldextest.Rdata", envir = D4STD0301)
load("/home/pvannameihome/PG0xxx03-01/dada2fxxx/xxx_allphyloseq.Rdata", envir = D4STD0301)

ls(D4STD0301, all.names = TRUE)

### load  ###############

### from from STDncodseq04-01
D4STD0401 <- new.env(parent = emptyenv())

### load up files
load(file = "/home/pvannameihome/PG0xxx04-01/dada2fxxx/xxx_phyloseq.Rdata", envir = D4STD0401)

ls(D4STD0401, all.names = TRUE)

testcolk = pal_simpsons("springfield", alpha = 0.75)(16)
show_col(testcolk)
show_col(testcolk[c(6, 2, 3, 1, 8,
                    11, 10, 16, 13, 12,
                    7, 4, 5, 9)])

wkcolpal = testcolk[c(1, 5, 6,
                      15, 2, 3, 9, 8,
                      11, 10, 16, 13, 12,
                      7, 4, 14)]

names(wkcolpal) = c("n0", "n4", "n6", 0:12)
show_col(wkcolpal)

### INPUT should be from
###  STDncodseq03  ###
##### set params
fpond = "pg"
ftiss = "16S-rmchloMTunc"

fminprop = "001"
fminoc =  "05"

#reset the number of OTUs after removed
nbrOTUs = dim(D4STD0301$fcount_tab)[1]

fileprefx = paste("202109-", fpond, "-",
                  ftiss, "-",
                  "Q1f230-ald01-minprop", fminprop, "minOc", fminoc, "OTUs", nbrOTUs, sep = "")

pg2.p1.1 = 1
pg2.p1.2 = 4
pg2.p1.3 = c(6, 7)

g1 = "g2" #pg
g2 = "p1" #pd
 
p.list.1 = list(pg2.p1.1, pg2.p1.1,pg2.p1.2)
p.list.2 = list(pg2.p1.2, pg2.p1.3, pg2.p1.3)

p.para.1 = c("1", "1", "4")
p.para.2 = c("4", "67", "67")

#### test for pg0102b1x
# pg0102b1x sigald_fw_vs_nw

### output directory within the working directory

for (i in 1:length(p.list.1)) {

f.grp = p.para.1[i]
n.grp = p.para.2[i]

para <- paste(fpond, "-ald01-pgg2-pdp1-w", f.grp, "-vs-w", n.grp, sep ="")
print(para)

outdir = paste("/home/pvannameihome/PG0102b1x-ald01/", para, sep ="")

if (file.exists(outdir)){
    setwd(file.path(outdir))
} else {
    dir.create(file.path(outdir))
    setwd(file.path(outdir))
}

phyq1PG0102b1x.pair <- paste("phyq1PG0102b1x.pgg2-pdp1w", f.grp, "-vs-w", n.grp, sep ="")

phyq1PG0102b1x.pair = subset_samples(D4STD0301$f2.PG0102b1x_noChloMTunc_physeq, pg == g1 & pd == g2 & wk %in% c(p.list.1[[i]], p.list.2[[i]]))

asvsums = taxa_sums(phyq1PG0102b1x.pair)

kepasv.v = names(asvsums[which(asvsums > 0)])

phyq1PG0102b1x.pair <- prune_taxa(kepasv.v, phyq1PG0102b1x.pair)

######  test aldex
samples = sample_data(phyq1PG0102b1x.pair)

samples$srtest = samples$wk  ##should be changed everytime by P Ay
samples$srtest = as.vector(samples$srtest)

samples$srtest[which(samples$wk %in% c(p.list.1[[i]]))] = "f"
samples$srtest[which(samples$wk %in% c(p.list.2[[i]]))] = "n"
samples$srtest = factor(samples$srtest)

sample_data(phyq1PG0102b1x.pair) = samples

f2.conds = sample_data(phyq1PG0102b1x.pair)[colnames(otu_table(phyq1PG0102b1x.pair))]$srtest

# estimate the distribution of CLR values
f2.x <- aldex.clr(as(otu_table(phyq1PG0102b1x.pair), "matrix"), f2.conds,
                  mc.samples=1000, verbose=FALSE)

# generate the expected CLR value for each OTU
# along with expected value of effect sizes
f2.e <- aldex.effect(f2.x, f2.conds,
                     include.sample.summary=TRUE,
                     verbose=FALSE, useMC = FALSE)

# calculate expected P values
f2.tP <- aldex.ttest(f2.x, f2.conds, paired.test=FALSE)

### f2relCZM.PG0102b1x_noChloMTunc_physeq.besthit.ldf
dim(D4STD0401$f2relCZM.PG0102b1x_noChloMTunc_physeq.besthit.ldf)
ori.besthit.ldf = D4STD0401$f2relCZM.PG0102b1x_noChloMTunc_physeq.besthit.ldf

phyq1PG0102b1x.pair.besthit.ldf = ori.besthit.ldf[which(ori.besthit.ldf$oID %in% kepasv.v), ]


f2.mapIDs = distinct(phyq1PG0102b1x.pair.besthit.ldf[c("OTUID", "nID", "oID")], OTUID, nID, oID, .keep_all= TRUE)

f2.x.ald_all <- data.frame(f2.tP, f2.e)

f2.x.ald_all$ffoldchanges = (1/2)^f2.x.ald_all$diff.btw
f2.x.ald_all$aID = row.names(f2.x.ald_all)

f2.x.ald_all.m = merge(f2.x.ald_all, f2.mapIDs, by.x = "aID", by.y = "oID", all.x = TRUE)
row.names(f2.x.ald_all.m) = f2.x.ald_all.m$aID
f2.x.ald_all = f2.x.ald_all.m

if (sum(f2.x.ald_all$we.eBH < 0.05 & f2.x.ald_all$wi.eBH < 0.05 &  f2.x.ald_all$effect < -1) == 0 & sum(f2.x.ald_all$we.eBH < 0.05 & f2.x.ald_all$wi.eBH < 0.05 &  f2.x.ald_all$effect > 1) == 0) {

save(phyq1PG0102b1x.pair,
     phyq1PG0102b1x.pair.besthit.ldf,     
     file = paste(fileprefx, "_pgg2-pdp1_sigald_fw", f.grp, "_vs_n", n.grp, ".Rdata", sep = ""))

} else {

# get high in f
f2.P.lowp.higheE = which(f2.x.ald_all$we.eBH < 0.05 & f2.x.ald_all$wi.eBH < 0.05 & f2.x.ald_all$effect < -1)
f2.P.lowp.higheE.OTUID = rownames(f2.x.ald_all[f2.P.lowp.higheE, ])
f2.x.ald_all.higheE = f2.x.ald_all[f2.P.lowp.higheE, ]


# get high in n
f2.P.lowp.higheN = which(f2.x.ald_all$we.eBH < 0.05 & f2.x.ald_all$wi.eBH < 0.05 & f2.x.ald_all$effect > 1)
f2.P.lowp.higheN.OTUID = rownames(f2.x.ald_all[f2.P.lowp.higheN, ])
f2.x.ald_all.higheN = f2.x.ald_all[f2.P.lowp.higheN, ]

tablef2.x.ald_all.higheE <- datatable(f2.x.ald_all.higheE, filter = 'top', options = list(paging = FALSE))
tablef2.x.ald_all.higheN <- datatable(f2.x.ald_all.higheN, filter = 'top', options = list(paging = FALSE))

saveWidget(tablef2.x.ald_all.higheE, file = paste(fileprefx, "-pgg2-pdp1-Dif-fw", f.grp, "_vs_nw", n.grp, "_FmN_abundance_table.html", sep = "")
)
saveWidget(tablef2.x.ald_all.higheN, file = paste(fileprefx, "-pgg2-pdp1-Dif-fw", f.grp, "_vs_nw", n.grp, "_NmF_abundance_table.html", sep = "")
)

save(phyq1PG0102b1x.pair,
     phyq1PG0102b1x.pair.besthit.ldf,
		 f2.x.ald_all,
		 f2.x.ald_all.higheE,
		 f2.x.ald_all.higheN,
		 f2.P.lowp.higheE,
		 f2.P.lowp.higheE.OTUID,
		 f2.P.lowp.higheN,
		 f2.P.lowp.higheN.OTUID,
		 file = paste(fileprefx, "-pgg2-pdp1_sigald_fw", f.grp, "_vs_nw", n.grp, ".Rdata", sep = ""))     

###############

### from from STDncodseq04-01

f2relCZM.PG0102b1x_noChloMTunc_physeq.besthit.ldf = ori.besthit.ldf
names(f2relCZM.PG0102b1x_noChloMTunc_physeq.besthit.ldf)
f2relCZM.PG0102b1x_noChloMTunc_physeq.besthit.ldf[, c("Sam_rep", "pg", "source", "wk", "tissue", "pd", "pdwk")]

#-----------
wklist = unique(f2relCZM.PG0102b1x_noChloMTunc_physeq.besthit.ldf$wk)

## check wk list order
wklist

# [1] 1  4  6  7  10 11
# Levels: 1 4 6 7 10 11

orwklist = wklist

f2relCZM.PG0102b1x_noChloMTunc_physeq.besthit.ldf$wk = factor(f2relCZM.PG0102b1x_noChloMTunc_physeq.besthit.ldf$wk, levels = orwklist)


f2relCZM.PG0102b1x_noChloMTunc_physeq.besthit.ldf$Sam_rep = factor(f2relCZM.PG0102b1x_noChloMTunc_physeq.besthit.ldf$Sam_rep)

#x.sampledat = sample_data(f2relCZM.PG03p3Q1_noChloMTunc_physeq.besthit)

pdfxwk.order = sort(unique(f2relCZM.PG0102b1x_noChloMTunc_physeq.besthit.ldf$pdwk))

by.pdfxwk = pdfxwk.order[c(1:5,
                           8:10, 
                           6:7)] #row.names(x.sampledat[order(x.sampledat$pdwk), ])
#check order of by.pdfxwk 
by.pdfxwk 
# > by.pdfxwk
# [1] "p1-1"  "p1-4"  "p1-6"  "p1-7"  "p2-1"  "p2-4"  "p2-6"  "p2-7"  "p2-10" "p2-11"


## forgraphic below

f2relCZM.PG0102b1x_noChloMTunc_physeq.besthit.ldf$pdwk = factor(f2relCZM.PG0102b1x_noChloMTunc_physeq.besthit.ldf$pdwk , levels = by.pdfxwk)

sam.pdwk = unique(f2relCZM.PG0102b1x_noChloMTunc_physeq.besthit.ldf[, c("Sam_rep","pdwk")])

sam.pdwk = sam.pdwk[
  with(sam.pdwk, order(pdwk, Sam_rep)),
  ]

samp.by.pdwk = as.vector(sam.pdwk$Sam_rep)

samp.by.pdwk = factor(samp.by.pdwk, levels = samp.by.pdwk)


ALLphyseq = D4STD0401$f2relCZM.PG0102b1x_noChloMTunc_physeq.besthit
nASVidrown = rownames(tax_table(ALLphyseq))

nASVidrown_oid = c()
for (iz in 1:(length(nASVidrown))){
  tr1 = ((strsplit(as.character(nASVidrown[iz]), ":", fixed = TRUE))[[1]])[1]
  nASVidrown_oid  = c(nASVidrown_oid, tr1)
  #c(nASVidrown_oid, ((strsplit(as.character(tr1), "-", fixed = TRUE))[[1]])[2])
}

tax_table(ALLphyseq) <- cbind(tax_table(ALLphyseq),
                              nASVidrown_oid)

colnames(tax_table(ALLphyseq)) <-
  c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "best_hit", "ASVID")

#### top ---

### lines to divide categories
colnames(otu_table(ALLphyseq))
sam.pdwk$pdwk
sam.pdwk$pdwk[1:4]  # p1-1
sam.pdwk$pdwk[5:8] # p1-4
sam.pdwk$pdwk[9:12] # p1-6
sam.pdwk$pdwk[13:16] # p1-7
sam.pdwk$pdwk[17:20] #  p2-1
sam.pdwk$pdwk[21:24] #  p2-4
sam.pdwk$pdwk[25:28] # p2-6
sam.pdwk$pdwk[29:32] # p2-7
sam.pdwk$pdwk[33:36] # p2-10
sam.pdwk$pdwk[37:40] # p2-11

xcapcat = c(4.5, 8.5, 
            12.5, 
            16.5,  ##p2
            20.5, 
            24.5,  
            28.5, 32.3,
            36.5)
					 
xcolxcapcat=  c("blue", "blue", "blue", 
                "orange", 
                "blue", "blue","blue",
                "blue", "blue")

## for ALL ASV in E
grpName_ALL= paste("-pgg2-pdp1-E_fw", f.grp, "_vs_nw", n.grp,"_ASV", sep = "")

grp_ALL <- f2.P.lowp.higheE.OTUID

SLXphyseq_ALL = phyloseq::subset_taxa(ALLphyseq, ASVID %in% grp_ALL)

titleASVID_ALL= paste(c(grpName_ALL), collapse = "\n")

colnames(sample_data(SLXphyseq_ALL))
sample_data(SLXphyseq_ALL) <- cbind(sample_data(SLXphyseq_ALL),
                                    (sample_data(SLXphyseq_ALL))[, "pd"])

samdatx = as.data.frame(sample_data(SLXphyseq_ALL))
samdatx$wk = factor(samdatx$wk, levels = orwklist)
sample_data(SLXphyseq_ALL) = samdatx


### print ballon plot

SLXphyseq.ALL.ldf = phy_to_ldf(SLXphyseq_ALL, transform.counts = NULL)

SLXphyseq.ALL.ldf$Abundance = SLXphyseq.ALL.ldf$Abundance * 100

SLXphyseq.ALL.ldf[which(SLXphyseq.ALL.ldf$Abundance < 1e-2), "Abundance"] = NA

SLXphyseq.ALL.ldf$pdfxwkSam_rep = paste(SLXphyseq.ALL.ldf$pd, SLXphyseq.ALL.ldf$wk,
                                        SLXphyseq.ALL.ldf$Sam_rep, sep = "-")				
																																		
lst.pdfxwkSam_rep  = sort(unique(SLXphyseq.ALL.ldf$pdfxwkSam_rep))

lst.pdfxwkSam_rep = lst.pdfxwkSam_rep[c(1:20, 29:40, 21:28)]

SLXphyseq.ALL.ldf$pdfxwkSam_rep = factor(SLXphyseq.ALL.ldf$pdfxwkSam_rep, levels = lst.pdfxwkSam_rep)

mycolW <- rgb(255, 255, 255, max = 255, alpha = 25, names = "W25")

SLXphyseq.ALL.ldf$taxclass = paste(SLXphyseq.ALL.ldf$Phylum,
                                   SLXphyseq.ALL.ldf$Class,  SLXphyseq.ALL.ldf$Order,
                                   SLXphyseq.ALL.ldf$Family,
                                   sep = "-")

ASVID.taxclass = unique(SLXphyseq.ALL.ldf[, c("ASVID","taxclass")])

ASVID.taxclass = ASVID.taxclass[
  with(ASVID.taxclass, order(taxclass, ASVID)),
  ]

ASVID.by.taxclass = as.vector(ASVID.taxclass$ASVID)

SLXphyseq.ALL.ldf$ASVID = factor(SLXphyseq.ALL.ldf$ASVID, levels = ASVID.by.taxclass)

gbp = ggballoonplot(SLXphyseq.ALL.ldf, y="ASVID", x= "pdfxwkSam_rep", size = "Abundance",
                    color = "taxclass", fill = mycolW,
                    title = titleASVID_ALL)

gbp = gbp + geom_vline(xintercept = xcapcat, linetype="dotted",
                       color =xcolxcapcat, size=1)

gbp= gbp+  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
gbp

ggsave(filename=paste(fileprefx, "-", grpName_ALL, "_balloons_yASV1.pdf", sep = ""), width=800, height=200, units="mm",
       limitsize = FALSE)

write.table(tax_table(SLXphyseq_ALL), file =paste(fileprefx, "-", grpName_ALL, "-top.tab", sep = ""), sep = "\t" )





## for ASV in "N"
## for ALL ASV in N
grpName_A= paste("-pgg2-pdp1-N_fw", f.grp, "_vs_nw", n.grp, "_ASV", sep = "")

grp_A <- f2.P.lowp.higheN.OTUID

SLXphyseq_A = phyloseq::subset_taxa(ALLphyseq, ASVID %in% grp_A)

titleASVID_A= paste(c(grpName_A), collapse = "\n")

colnames(sample_data(SLXphyseq_A))
sample_data(SLXphyseq_A) <- cbind(sample_data(SLXphyseq_A),
                                    (sample_data(SLXphyseq_A))[, "pd"])

samdatx = as.data.frame(sample_data(SLXphyseq_A))
samdatx$wk = factor(samdatx$wk, levels = orwklist)
sample_data(SLXphyseq_A) = samdatx


### print ballon plot

SLXphyseq.A.ldf = phy_to_ldf(SLXphyseq_A, transform.counts = NULL)

SLXphyseq.A.ldf$Abundance = SLXphyseq.A.ldf$Abundance * 100

SLXphyseq.A.ldf[which(SLXphyseq.A.ldf$Abundance < 1e-2), "Abundance"] = NA

SLXphyseq.A.ldf$pdfxwkSam_rep = paste(SLXphyseq.A.ldf$pd, SLXphyseq.A.ldf$wk,
                                  SLXphyseq.A.ldf$Sam_rep, sep = "-")
																	
lst.pdfxwkSam_rep = sort(unique(SLXphyseq.A.ldf$pdfxwkSam_rep))

lst.pdfxwkSam_rep = lst.pdfxwkSam_rep[c(1:20, 29:40, 21:28)]

SLXphyseq.A.ldf$pdfxwkSam_rep = factor(SLXphyseq.A.ldf$pdfxwkSam_rep, levels = lst.pdfxwkSam_rep)

mycolW <- rgb(255, 255, 255, max = 255, alpha = 25, names = "W25")

SLXphyseq.A.ldf$taxclass = paste(SLXphyseq.A.ldf$Phylum,
                                   SLXphyseq.A.ldf$Class,
                                   SLXphyseq.A.ldf$Order,
                                   SLXphyseq.A.ldf$Family,
                                   sep = "-")

ASVID.taxclass = unique(SLXphyseq.A.ldf[, c("ASVID","taxclass")])

ASVID.taxclass = ASVID.taxclass[
  with(ASVID.taxclass, order(taxclass, ASVID)),
  ]

ASVID.by.taxclass = as.vector(ASVID.taxclass$ASVID)

SLXphyseq.A.ldf$ASVID = factor(SLXphyseq.A.ldf$ASVID, levels = ASVID.by.taxclass)

gbp = ggballoonplot(SLXphyseq.A.ldf, y="ASVID", x= "pdfxwkSam_rep", size = "Abundance",
              color = "taxclass", fill = mycolW,
              title = titleASVID_A)

gbp = gbp + geom_vline(xintercept = xcapcat, linetype="dotted",
                       color =xcolxcapcat, size=1)
	
gbp= gbp+  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
gbp

ggsave(filename=paste(fileprefx, "-", grpName_A, "_balloons_yASV1.pdf", sep = ""), width=880, height=200, units="mm",
       limitsize = FALSE)

write.table(tax_table(SLXphyseq_A), file =paste(fileprefx, "-", grpName_A, "-top.tab", sep = ""), sep = "\t" )

}

}

