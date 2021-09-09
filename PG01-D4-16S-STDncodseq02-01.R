
# commands to load libraries and datasets
library(phyloseq)
library(qiime2R)
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
library(ggsci)
library(scales)
library(ggplot2)
library(DESeq2)
library(dendextend)
library(microbiome) # data analysis and visualisation
library(microbiomeutilities) # some utility tools
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(DT) # interactive tables in html and markdown
library(data.table) # alternative to data.frame

##### STDnCoda_pipeline

### READING IN OUR DATA  ## need to set this
setwd("/home/sylvester/test/PG01-D4-16S-STDncodseq02-01")

##### set params
fpond = "pg01-D4"
ftiss = "16S-rmchloMTunc"
prop = "Q1f230-STDnCoda01-minprop"
prop1 = "Q1f230-STDnCoda01-cutoffs"

input_qza = "dada2f230_pg01ds0123416s_tab.qza"
input_tab = "dada2f230_pg01234_rep_silvaV3V4-97_tax.tab"
input_chlo = "dada2f230_pg01234_rep_silvaV3V4-97_tax.Chlo"
input_mit = "dada2f230_pg01234_rep_silvaV3V4-97_tax.mitch"
input_cul = "dada2f230_pg01234_rep_silvaV3V4-97_tax.ricket.uncul"
input_tsv = "pg01-ds01234-metadata.tsv"


# This command should work as-is with Greengenes.
phyq1PG01 <- qza_to_phyloseq(features=input_qza, metadata=input_tsv )

count_tab <- data.frame(otu_table(phyq1PG01))

tax_tab <- as.matrix(read.table(input_tab, header=T, row.names=1, check.names=F, na.strings="", sep="\t"))

tax_tab.chlo <- as.matrix(read.table(input_chlo, header=F, row.names=1, check.names=F, na.strings="", sep="\t"))

tax_tab.chlo.asvid = row.names(tax_tab.chlo)


### mitochondria related
tax_tab.MT <- as.matrix(read.table(input_mit, header=F, row.names=1, check.names=F, na.strings="", sep="\t"))

tax_tab.MT.asvid = row.names(tax_tab.MT)

#tax_tab.rick.unc <- as.matrix(read.table(input_cul, header=F, row.names=1, check.names=F, na.strings="", sep="\t"))

#tax_tab.rick.unc.asvid = row.names(tax_tab.rick.unc)

#tax_tab.MTunc.asvid = c(tax_tab.MT.asvid, tax_tab.rick.unc.asvid)

tax_tab.MTunc.asvid = c(tax_tab.MT.asvid)

sample_info_tab = sample_data(phyq1PG01)


### output directory within the working directory
outdir = "/home/sylvester/test/PG01-D4-16S-STDncodseq02-01/dada2f230/"

if (file.exists(outdir)){
    setwd(file.path(outdir))
} else {
    dir.create(file.path(outdir))
    setwd(file.path(outdir))
}

testcolk = pal_simpsons("springfield", alpha = 0.75)(16)

wkcolpal = testcolk[c(15, 2, 3, 9, 8, 11, 10, 16, 13, 12, 7, 4, 14)]

names(wkcolpal) = 0:12


## remove chlo
count_tab.chlo.rowid = which(row.names(count_tab) %in% tax_tab.chlo.asvid)
filt_count_tab.nochlo = count_tab[-count_tab.chlo.rowid , ]
tax_tab.nochlo = tax_tab[-which(row.names(tax_tab) %in% tax_tab.chlo.asvid), ]

## check Chloroplast
which(tax_tab.nochlo[, 3] == "Chloroplast")

## remove MTunc
count_tab.MTunc.rowid = which(row.names(filt_count_tab.nochlo) %in% tax_tab.MTunc.asvid)
filt_count_tab = filt_count_tab.nochlo[-count_tab.MTunc.rowid , ]
tax_tab.nochloMTunc = tax_tab.nochlo[-which(row.names(tax_tab.nochlo) %in% tax_tab.MTunc.asvid), ]

#check Mitochondria
which(tax_tab.nochloMTunc[, 3] == "Mitochondria")

filt_sample_info_tab = sample_info_tab

### get column_sum

totalreads = count_tab %>%
  replace(is.na(.), 0) %>%
  summarise_all(list(sum))

totalreads.nochlo = filt_count_tab.nochlo %>%
  replace(is.na(.), 0) %>%
  summarise_all(list(sum))

totalreads.nochloMTunc = filt_count_tab %>%
  replace(is.na(.), 0) %>%
  summarise_all(list(sum))

totalreads.ofchloMTunc = totalreads - totalreads.nochloMTunc

totalreads.ofchlo = totalreads - totalreads.nochlo

totalreads.ofMTunc = totalreads.ofchloMTunc - totalreads.ofchlo

## count ASVs
## head(tax_tab)
## head(tax_tab.chlo)
totalASV = dim(tax_tab)[1]
totalASV.chlo = dim(tax_tab.chlo)[1]
totalASV.MTunc = length(tax_tab.MTunc.asvid)

### get min occur
2/dim(filt_sample_info_tab)[1]

#[1] 0.02531646 for fminoc

############################## filter for frequency   filter-02
## need to set min.occurrence=
## and  min.prop=


fcount_tab <- codaSeq.filter(count_tab, min.reads=5000, min.prop=0.01, min.occurrence=0.02, samples.by.row=FALSE)

## for ALL using this.
fminprop = "01"
fminoc =  "02"

##nbrtotal OTU
nbrOTUs = dim(fcount_tab)[1]    #255

### remove chlo
fcount_tab.chlo.rowid = which(row.names(fcount_tab) %in% tax_tab.chlo.asvid)
filt_fcount_tab = fcount_tab[-fcount_tab.chlo.rowid , ]
chlor_fcount_tab = fcount_tab[fcount_tab.chlo.rowid , ]

fcount_tab.chloMTunc.rowid = which(row.names(filt_fcount_tab) %in% tax_tab.MTunc.asvid)

if(length(fcount_tab.chloMTunc.rowid) > 0){
  filt_fcount_tab = filt_fcount_tab[-fcount_tab.chloMTunc.rowid , ]
}else{
  filt_fcount_tab = filt_fcount_tab
}

MTunc_fcount_tab = filt_fcount_tab[fcount_tab.chloMTunc.rowid , ]


##count ASVs
f.totalASV = dim(fcount_tab)[1]
f.totalASV.chlo = dim(chlor_fcount_tab)[1]
f.totalASV.MTunc= dim(MTunc_fcount_tab)[1]

f_totalreads = fcount_tab %>%
  replace(is.na(.), 0) %>%
  summarise_all(list(sum))

f_totalreads.nochloMTunc  = filt_fcount_tab %>%
  replace(is.na(.), 0) %>%
  summarise_all(list(sum))

f_totalreads.ofchlo  = chlor_fcount_tab %>%
  replace(is.na(.), 0) %>%
  summarise_all(list(sum))

f_totalreads.ofMTunc  = MTunc_fcount_tab %>%
  replace(is.na(.), 0) %>%
  summarise_all(list(sum))


####
fcount_tab <- codaSeq.filter(count_tab, min.reads=5000,
                             min.prop=0.001, min.occurrence=0.02, samples.by.row=FALSE)
## for ALL using this.
fminprop = "001"

##nbrtotal OTU
nbrOTUs = dim(fcount_tab)[1]    #1195
#######################################

###########
## file prefix
fileprefx = paste("202102-", fpond, "-",
                  ftiss, "-",
                  prop, "-", fminprop, "minOc", fminoc, "OTUs", nbrOTUs, sep = "")

fileprefx
###########################################

### remove chlo
fcount_tab.chlo.rowid = which(row.names(fcount_tab) %in% tax_tab.chlo.asvid)
filt_fcount_tab = fcount_tab[-fcount_tab.chlo.rowid , ]
chlor_fcount_tab = fcount_tab[fcount_tab.chlo.rowid , ]

fcount_tab.chloMTunc.rowid = which(row.names(filt_fcount_tab) %in% tax_tab.MTunc.asvid)

MTunc_fcount_tab = filt_fcount_tab[fcount_tab.chloMTunc.rowid , ]

if(length(fcount_tab.chloMTunc.rowid) > 0){
  filt_fcount_tab = filt_fcount_tab[-fcount_tab.chloMTunc.rowid , ]
}else{
  filt_fcount_tab = filt_fcount_tab
}


##count ASVs
f01.totalASV = dim(fcount_tab)[1]
f01.totalASV.chlo = dim(chlor_fcount_tab)[1]
f01.totalASV.MTunc = dim(MTunc_fcount_tab)[1]

### get column_sum

f01_totalreads = fcount_tab %>%
  replace(is.na(.), 0) %>%
  summarise_all(list(sum))

f01_totalreads.nochloMTunc  = filt_fcount_tab %>%
  replace(is.na(.), 0) %>%
  summarise_all(list(sum))

f01_totalreads.ofchlo  = chlor_fcount_tab %>%
  replace(is.na(.), 0) %>%
  summarise_all(list(sum))

f01_totalreads.ofMTunc  = MTunc_fcount_tab %>%
  replace(is.na(.), 0) %>%
  summarise_all(list(sum))

#### p.0001
fcount_tab <- codaSeq.filter(count_tab, min.reads=5000,
                             min.prop=0.0001, min.occurrence=0.02, samples.by.row=FALSE)
## for ALL using this.
fminprop = "0001"

##nbrtotal OTU
nbrOTUs = dim(fcount_tab)[1]    #2699


### remove chlo
fcount_tab.chlo.rowid = which(row.names(fcount_tab) %in% tax_tab.chlo.asvid)
filt_fcount_tab = fcount_tab[-fcount_tab.chlo.rowid , ]
chlor_fcount_tab = fcount_tab[fcount_tab.chlo.rowid , ]

fcount_tab.chloMTunc.rowid = which(row.names(filt_fcount_tab) %in% tax_tab.MTunc.asvid)
MTunc_fcount_tab = filt_fcount_tab[fcount_tab.chloMTunc.rowid , ]

if(length(fcount_tab.chloMTunc.rowid) > 0){
  filt_fcount_tab = filt_fcount_tab[-fcount_tab.chloMTunc.rowid , ]
}else{
  filt_fcount_tab = filt_fcount_tab
}

##count ASVs
f001.totalASV = dim(fcount_tab)[1]
f001.totalASV.chlo = dim(chlor_fcount_tab)[1]
f001.totalASV.MTunc = dim(MTunc_fcount_tab)[1]
### get column_sum

f001_totalreads = fcount_tab %>%
  replace(is.na(.), 0) %>%
  summarise_all(list(sum))

f001_totalreads.nochloMTunc  = filt_fcount_tab %>%
  replace(is.na(.), 0) %>%
  summarise_all(list(sum))

f001_totalreads.ofchlo  = chlor_fcount_tab %>%
  replace(is.na(.), 0) %>%
  summarise_all(list(sum))

f001_totalreads.ofMTunc  = MTunc_fcount_tab %>%
  replace(is.na(.), 0) %>%
  summarise_all(list(sum))

cmb.ch = data.frame(all.totalr = as(t(totalreads[1, ]), "vector"),
                    all.totalr.nochMTunc = as(t(totalreads.nochloMTunc[1, ]), "vector"),
                    all.totalr.ofch = as(t(totalreads.ofchlo[1, ]), "vector"),
                    all.totalr.ofMTunc = as(t(totalreads.ofMTunc[1, ]), "vector"),

                    #           all.p.ofch = as(t(totalreads.ofchlo/totalreads,
                    p1.totalr = as(t(f_totalreads[1, ]), "vector"),
                    p1.totalr.nochMTunc = as(t(f_totalreads.nochloMTunc[1, ]), "vector"),
                    p1.totalr.ofch = as(t(f_totalreads.ofchlo[1, ]), "vector"),
                    p1.totalr.ofMTunc = as(t(f_totalreads.ofMTunc[1, ]), "vector"),

                    #           p1.p.ofch = f_totalreads.ofchlo/f_totalreads,
                    p01.totalr = as(t(f01_totalreads[1, ]), "vector"),
                    p01.totalr.nochMTunc = as(t(f01_totalreads.nochloMTunc[1, ]), "vector"),
                    p01.totalr.ofch = as(t(f01_totalreads.ofchlo[1, ]), "vector"),
                    p01.totalr.ofMTunc = as(t(f01_totalreads.ofMTunc[1, ]), "vector"),

                    p001.totalr = as(t(f001_totalreads[1, ]), "vector"),
                    p001.totalr.nochMTunc = as(t(f001_totalreads.nochloMTunc[1, ]), "vector"),
                    p001.totalr.ofch = as(t(f001_totalreads.ofchlo[1, ]), "vector"),
                    p001.totalr.ofMTunc = as(t(f001_totalreads.ofMTunc[1, ]), "vector")

                    #           p01.p.ofch = f01_totalreads.ofchlo/f01_totalreads
)

row.names(cmb.ch) = names(totalreads)


### plot percent chloroplast
p.cmb.ch = data.frame(all = cmb.ch$all.totalr.ofch/cmb.ch$all.totalr,
                      p1 = cmb.ch$p1.totalr.ofch/cmb.ch$p1.totalr,
                      p01 =cmb.ch$p01.totalr.ofch/cmb.ch$p01.totalr,
                      p001 =cmb.ch$p001.totalr.ofch/cmb.ch$p001.totalr)

row.names(p.cmb.ch) = row.names(cmb.ch)
p.cmb.ch$sampid = row.names(p.cmb.ch)
p.cmb.ch.lng = gather(p.cmb.ch, cutoff, pchlor, all:p001, factor_key=TRUE)

df.sample_info_tab = as(sample_info_tab, "data.frame")
df.sample_info_tab$sampid = row.names(df.sample_info_tab)

p.cmb.ch.saminfo = merge(p.cmb.ch.lng, df.sample_info_tab,
                         by.x = "sampid", by.y = "sampid" )

p.cmb.ch.saminfo$idno = paste(p.cmb.ch.saminfo$shrimpid,
                              p.cmb.ch.saminfo$site, p.cmb.ch.saminfo$filter,
                              sep="-")

p.cmb.ch.saminfo$cutoff = factor(p.cmb.ch.saminfo$cutoff,levels(p.cmb.ch.saminfo$cutoff)[
  c(1,4,3,2)])


psamp = ggscatter(data = p.cmb.ch.saminfo, x="wk", y="pchlor",
                  #                  xlab=xlab, ylab = ylab,
                  size = 3,
                  facet.by = c("source", "cutoff"),
                  #                  main = paste(fileprefx, "PCA samples", sep = "\n"),
                  #                  xlim = c(min(muldifsample$pc1) -5, max(muldifsample$pc1) +5    ),
                  #                  ylim = c(min(muldifsample$pc2) -5, max(muldifsample$pc2) +5    ),
                  color="wk",
                  shape = "source",
                  position = position_jitter(width = 0.2, height = 0),
                  label = "idno",
                  repel = TRUE ,
                  #                  label.select = c("f8", "f6","n5"),
                  #                  ellipse = TRUE,
                  palette = wkcolpal
)

psamp +
  #  stat_conf_ellipse(aes(color = shrimp, fill = pg01), alpha = 0.1, geom = "polygon") +
  #  geom_jitter(width = 0.2)+
  facet_grid(rows = vars(cutoff),cols = vars(source), scales="free") +
  #  scale_y_continuous(trans="log10") +
  geom_hline(yintercept=0, linetype=2, size=1, color=rgb(0,0,0,0.3)) +
  geom_vline(xintercept=0, linetype=2, size=1, color=rgb(0,0,0,0.3))


ggsave(filename=paste("202102-", fpond, "-",
                      ftiss, "-",
                      prop1, "_pctchloro_label.pdf", sep = ""),
                                                                width=390, height=210, units="mm")

psamp = ggscatter(data = p.cmb.ch.saminfo, x="wk", y="pchlor",
                  #                  xlab=xlab, ylab = ylab,
                  size = 3,
                  facet.by = c("source", "cutoff"),
                  #                  main = paste(fileprefx, "PCA samples", sep = "\n"),
                  #                  xlim = c(min(muldifsample$pc1) -5, max(muldifsample$pc1) +5    ),
                  #                  ylim = c(min(muldifsample$pc2) -5, max(muldifsample$pc2) +5    ),
                  color="wk",
                  shape = "source",
                  position = position_jitter(width = 0.2, height = 0),
                  #label = "idno",
                  repel = TRUE ,
                  #                  label.select = c("f8", "f6","n5"),
                  #                  ellipse = TRUE,
                  palette = wkcolpal
)

psamp +
  #  stat_conf_ellipse(aes(color = shrimp, fill = pg01), alpha = 0.1, geom = "polygon") +
  #  geom_jitter(width = 0.2)+
  facet_grid(rows = vars(cutoff),cols = vars(source), scales="free") +
  #  scale_y_continuous(trans="log10") +
  geom_hline(yintercept=0, linetype=2, size=1, color=rgb(0,0,0,0.3)) +
  geom_vline(xintercept=0, linetype=2, size=1, color=rgb(0,0,0,0.3))

ggsave(filename=paste("202102-", fpond, "-",
                      ftiss, "-",
                      prop1, "_pctchloro_xbel.pdf", sep = ""),
                                                                width=390, height=210, units="mm")

### plot percent MTunc
p.cmb.MTunc = data.frame(all = cmb.ch$all.totalr.ofMTunc/cmb.ch$all.totalr,
                      p1 = cmb.ch$p1.totalr.ofMTunc/cmb.ch$p1.totalr,
                      p01 =cmb.ch$p01.totalr.ofMTunc/cmb.ch$p01.totalr,
                      p001 =cmb.ch$p001.totalr.ofMTunc/cmb.ch$p001.totalr)

row.names(p.cmb.MTunc) = row.names(cmb.ch)
p.cmb.MTunc$sampid = row.names(p.cmb.MTunc)
p.cmb.MTunc.lng = gather(p.cmb.MTunc, cutoff, pMTunc, all:p001, factor_key=TRUE)

df.sample_info_tab = as(sample_info_tab, "data.frame")
df.sample_info_tab$sampid = row.names(df.sample_info_tab)

p.cmb.MTunc.saminfo = merge(p.cmb.MTunc.lng, df.sample_info_tab,
                                        by.x = "sampid", by.y = "sampid" )

p.cmb.MTunc.saminfo$idno = paste(p.cmb.MTunc.saminfo$shrimpid,
                                 p.cmb.MTunc.saminfo$site, p.cmb.MTunc.saminfo$filter,
                                 sep="-")

p.cmb.MTunc.saminfo$cutoff = factor(p.cmb.MTunc.saminfo$cutoff,levels(p.cmb.MTunc.saminfo$cutoff)[c(1,4,3,2)])

psamp = ggscatter(data = p.cmb.MTunc.saminfo, x="wk", y="pMTunc",
                  #                  xlab=xlab, ylab = ylab,
                  size = 3,
                  facet.by = c("source", "cutoff"),
                  #                  main = paste(fileprefx, "PCA samples", sep = "\n"),
                  #                  xlim = c(min(muldifsample$pc1) -5, max(muldifsample$pc1) +5    ),
                  #                  ylim = c(min(muldifsample$pc2) -5, max(muldifsample$pc2) +5    ),
                  color="wk",
                  shape = "source",
                  position = position_jitter(width = 0.2, height = 0),
                  label = "idno",
                  repel = TRUE ,
                  #                  label.select = c("f8", "f6","n5"),
                  #                  ellipse = TRUE,
                  palette = wkcolpal
)

psamp +
  #  stat_conf_ellipse(aes(color = shrimp, fill = pg01), alpha = 0.1, geom = "polygon") +
  #  geom_jitter(width = 0.2)+
  facet_grid(rows = vars(cutoff),cols = vars(source), scales="free") +
  #  scale_y_continuous(trans="log10") +
  geom_hline(yintercept=0, linetype=2, size=1, color=rgb(0,0,0,0.3)) +
  geom_vline(xintercept=0, linetype=2, size=1, color=rgb(0,0,0,0.3))

ggsave(filename=paste("202102-", fpond, "-",
                      ftiss, "-",
                      prop1, "_pctMTunc_label.pdf", sep = ""),
                                                          width=390, height=210, units="mm")

psamp = ggscatter(data = p.cmb.MTunc.saminfo, x="wk", y="pMTunc",
                  #                  xlab=xlab, ylab = ylab,
                  size = 3,
                  facet.by = c("source", "cutoff"),
                  #                  main = paste(fileprefx, "PCA samples", sep = "\n"),
                  #                  xlim = c(min(muldifsample$pc1) -5, max(muldifsample$pc1) +5    ),
                  #                  ylim = c(min(muldifsample$pc2) -5, max(muldifsample$pc2) +5    ),
                  color="wk",
                  shape = "source",
                  position = position_jitter(width = 0.2, height = 0),
                  #label = "idno",
                  repel = TRUE ,
                  #                  label.select = c("f8", "f6","n5"),
                  #                  ellipse = TRUE,
                  palette = wkcolpal
)

psamp +
  #  stat_conf_ellipse(aes(color = shrimp, fill = pg01), alpha = 0.1, geom = "polygon") +
  #  geom_jitter(width = 0.2)+
  facet_grid(rows = vars(cutoff),cols = vars(source), scales="free") +
  #  scale_y_continuous(trans="log10") +
  geom_hline(yintercept=0, linetype=2, size=1, color=rgb(0,0,0,0.3)) +
  geom_vline(xintercept=0, linetype=2, size=1, color=rgb(0,0,0,0.3))

ggsave(filename=paste("202102-", fpond, "-",
                      ftiss, "-",
                      prop1, "_pctMTunc_xbel.pdf", sep = ""),
                                                          width=390, height=210, units="mm")

##   writeout info
write.table(p.cmb.ch.saminfo, file=paste("202102-", fpond, "-",
                                          ftiss, "-",
                                          prop1, "_pctchloro.tab", sep = ""),
                                                                                    sep = "\t", quote = FALSE)

write.table(p.cmb.MTunc.saminfo, file=paste("202102-", fpond, "-",
                                         ftiss, "-",
                                         prop1, "_pctMTunc.tab", sep = ""),
                                                                                     sep = "\t", quote = FALSE)

write.table(cmb.ch, file=paste("202102-", fpond, "-",
                                         ftiss, "-",
                                         prop1, "_ReadchloroMTunc.tab", sep = ""),
                                                                                     sep = "\t", quote = FALSE)

### count ASVs
cmb.ASVch = data.frame(ta = c(totalASV, f.totalASV, f01.totalASV,  f001.totalASV  ),
                       tch = c(totalASV.chlo, f.totalASV.chlo, f01.totalASV.chlo,  f001.totalASV.chlo),
                       tMTunc = c(totalASV.MTunc, f.totalASV.MTunc, f01.totalASV.MTunc,  f001.totalASV.MTunc)
)

cmb.ASVch$tnch = cmb.ASVch$ta - cmb.ASVch$tch
cmb.ASVch$tnchMTunc = cmb.ASVch$ta - cmb.ASVch$tch - cmb.ASVch$tMTunc
cmb.ASVch$pcnt = c("all", "1%", "0.1%", "0.01%")
row.names(cmb.ASVch) = c("all", "p1", "p01", "p001")

write.table(cmb.ASVch, file=paste("202102-", fpond, "-",
                               ftiss, "-",
                               prop1, "_ASVchloroMTunc.tab", sep = ""),
                                                                                 sep = "\t", quote = FALSE)


#### P001
##  "0.1%"
####
### need to change min.prop=0.001, min.occurrence  based on the above

fcount_tab <- codaSeq.filter(count_tab, min.reads=5000,
                             min.prop=0.001, min.occurrence=0.02, samples.by.row=FALSE)
## for ALL using this.
fminprop = "001"

##nbrtotal OTU
nbrOTUs = dim(fcount_tab)[1]    #1195

### remove chlo
fcount_tab.chlo.rowid = which(row.names(fcount_tab) %in% tax_tab.chlo.asvid)
filt_fcount_tab = fcount_tab[-fcount_tab.chlo.rowid , ]
chlor_fcount_tab = fcount_tab[fcount_tab.chlo.rowid , ]

fcount_tab.chloMTunc.rowid = which(row.names(filt_fcount_tab) %in% tax_tab.MTunc.asvid)
MTunc_fcount_tab = filt_fcount_tab[fcount_tab.chloMTunc.rowid , ]
if(length(fcount_tab.chloMTunc.rowid) > 0){
  filt_fcount_tab = filt_fcount_tab[-fcount_tab.chloMTunc.rowid , ]
}else{
  filt_fcount_tab = filt_fcount_tab
}

filt_tax_tab = tax_tab[row.names(filt_fcount_tab), ]

chlor_tax_tab = tax_tab[row.names(chlor_fcount_tab), ]

MTunc_tax_tab = tax_tab[row.names(MTunc_fcount_tab), ]

#### phyloseq


# replace 0 values with an estimate  ######### f.n0.count_tab
f.n0.count_tab <- cmultRepl(t(fcount_tab), method="CZM", label=0)

#output message
#No. corrected values:  27

# generate the CLR values for plotting later
f.n0.clr <- codaSeq.clr(f.n0.count_tab)

####Returns a matrix of center log-ratio transformed data with samples by row.
##Each value is equivalent to log(x/gx) where gx is the geometic mean of the row vector X.


###### reconstruct phyloseq again with the raw count and CZM method (replace 0 with an estimated).
#Transform into matrixes otu and tax tables
otu_mat <- as.matrix(fcount_tab)
tax_mat <- tax_tab[row.names(fcount_tab),]

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(filt_sample_info_tab)  ### change to filter for ponds/tissues

### please also change the prefix file e.g., PG01p301 to pg0104f230Q1 for all the corresponding objects / variables

f2.PG01p3Q1physeq <- phyloseq(OTU, TAX, samples)
### phyloseq with the counts.

f2relCZM.PG01p3Q1physeq = f2.PG01p3Q1physeq

otu_mat <- as.matrix(t(f.n0.count_tab))  ### transposed f.n0

otu_table(f2relCZM.PG01p3Q1physeq) = otu_table(otu_mat, taxa_are_rows = TRUE)

f2relCZM.PG01p3Q1physeq ## phyloseq with estimated 0 count from CZM method

### codaSeq.clr results
f2relCodaCRL.PG01p3Q1physeq = f2.PG01p3Q1physeq

otu_mat <- as.matrix(t(f.n0.clr))  ### transposed f.n0.clr

otu_table(f2relCodaCRL.PG01p3Q1physeq) = otu_table(otu_mat, taxa_are_rows = TRUE)

## phyloseq with CodaCRL


### microbiome transform cleer
f2rel.microCLR.PG01p3Q1physeq = microbiome::transform(f2.PG01p3Q1physeq, 'clr')

#### remove chlo ONLY
f2count_tab.chlo.rowid = which(row.names(otu_table(f2.PG01p3Q1physeq)) %in%
                                 c(tax_tab.chlo.asvid ))

chlor_otu_tab = subset(otu_table(f2.PG01p3Q1physeq),
                       rownames(otu_table(f2.PG01p3Q1physeq)) %in%
                       c(tax_tab.chlo.asvid ))

f2.PG01p3Q1_chlorONLY_physeq = merge_phyloseq(chlor_otu_tab,
                                          tax_table(f2.PG01p3Q1physeq),
                                          sample_data(f2.PG01p3Q1physeq))

#### get CMZ

chlor_otu_tab = subset(otu_table(f2relCZM.PG01p3Q1physeq),
                       rownames(otu_table(f2relCZM.PG01p3Q1physeq)) %in%
                       c(tax_tab.chlo.asvid))

f2relCZM.PG01p3Q1_chlorONLY_physeq = merge_phyloseq(chlor_otu_tab,
                                                        tax_table(f2relCZM.PG01p3Q1physeq),
                                                        sample_data(f2relCZM.PG01p3Q1physeq))


## print tax
plot_bar(f2relCZM.PG01p3Q1_chlorONLY_physeq, fill="Family") + theme_bw() +  theme(axis.text.x = element_text(angle = 90))

ggsave(filename=paste(fileprefx, "-chloONLY", "_familyAbundance.pdf", sep = ""), width=420, height=230, units="mm")

#### remove MTunc ONLY
f2count_tab.MTunc.rowid = which(row.names(otu_table(f2.PG01p3Q1physeq)) %in% c( tax_tab.MTunc.asvid ))

MTunc_otu_tab = subset(otu_table(f2.PG01p3Q1physeq),
                       rownames(otu_table(f2.PG01p3Q1physeq)) %in%
                       c(tax_tab.MTunc.asvid ))

f2.PG01p3Q1_MTunc_physeq = merge_phyloseq(MTunc_otu_tab,
                                          tax_table(f2.PG01p3Q1physeq),
                                          sample_data(f2.PG01p3Q1physeq))


#### get CMZ

MTunc_otu_tab = subset(otu_table(f2relCZM.PG01p3Q1physeq),
                       rownames(otu_table(f2relCZM.PG01p3Q1physeq)) %in%
                       c(tax_tab.MTunc.asvid ))

f2relCZM.PG01p3Q1_MTuncONLY_physeq = merge_phyloseq(MTunc_otu_tab,
                                                tax_table(f2relCZM.PG01p3Q1physeq),
                                                sample_data(f2relCZM.PG01p3Q1physeq))


## print tax
plot_bar(f2relCZM.PG01p3Q1_MTuncONLY_physeq, fill="Genus") +theme_bw() +  theme(axis.text.x = element_text(angle = 90))

ggsave(filename=paste(fileprefx, "-MTuncONLY", "_GenusAbundance.pdf", sep = ""), width=420, height=230, units="mm")

### remove chlo & MTunc
f2count_tab.chlo.rowid = which(row.names(otu_table(f2.PG01p3Q1physeq)) %in% c(tax_tab.chlo.asvid, tax_tab.MTunc.asvid ))

chlor_otu_tab = subset(otu_table(f2.PG01p3Q1physeq),
                       rownames(otu_table(f2.PG01p3Q1physeq)) %in%
                       c(tax_tab.chlo.asvid, tax_tab.MTunc.asvid ))

f2.PG01p3Q1_chlor_physeq = merge_phyloseq(chlor_otu_tab,
                                              tax_table(f2.PG01p3Q1physeq),
                                              sample_data(f2.PG01p3Q1physeq))

filt_otu_tab = subset(otu_table(f2.PG01p3Q1physeq),
                      !(rownames(otu_table(f2.PG01p3Q1physeq)) %in%
                          c(tax_tab.chlo.asvid, tax_tab.MTunc.asvid )))

f2.PG01p3Q1_NonChl_physeq = merge_phyloseq(filt_otu_tab,
                                               tax_table(f2.PG01p3Q1physeq),
                                               sample_data(f2.PG01p3Q1physeq))


#### get CMZ f2relCZM.PG01p3Q1physeq

chlor_otu_tab = subset(otu_table(f2relCZM.PG01p3Q1physeq),
                       rownames(otu_table(f2relCZM.PG01p3Q1physeq)) %in%
                       c(tax_tab.chlo.asvid, tax_tab.MTunc.asvid ))
f2relCZM.PG01p3Q1_chlor_physeq = merge_phyloseq(chlor_otu_tab,
                                                    tax_table(f2relCZM.PG01p3Q1physeq),
                                                    sample_data(f2relCZM.PG01p3Q1physeq))

filt_otu_tab = subset(otu_table(f2relCZM.PG01p3Q1physeq),
                      !(rownames(otu_table(f2relCZM.PG01p3Q1physeq)) %in%
                          c(tax_tab.chlo.asvid, tax_tab.MTunc.asvid )))

f2relCZM.PG01p3Q1_NonChl_physeq = merge_phyloseq(filt_otu_tab,
                                                     tax_table(f2relCZM.PG01p3Q1physeq),
                                                     sample_data(f2relCZM.PG01p3Q1physeq))

#### save dataset

##### need to save objects
#### save phyloseq objects
save(
  fcount_tab,
  f2.PG01p3Q1physeq,
  f2rel.microCLR.PG01p3Q1physeq,
  f2relCodaCRL.PG01p3Q1physeq,
  f2relCZM.PG01p3Q1physeq,
  filt_count_tab,
  filt_sample_info_tab,
  phyq1PG01,
  tax_tab,
  tax_tab.chlo,
  tax_tab.nochlo,
  tax_tab.MT,
  #tax_tab.rick.unc,

  f2.PG01p3Q1_chlorONLY_physeq,
  f2.PG01p3Q1_MTunc_physeq,
  f2.PG01p3Q1_NonChl_physeq,
  f2relCZM.PG01p3Q1_chlorONLY_physeq,
  f2relCZM.PG01p3Q1_MTuncONLY_physeq,
  f2relCZM.PG01p3Q1_NonChl_physeq,
  tax_tab.chlo.asvid,
  tax_tab.MTunc.asvid,
  tax_tab.MT.asvid,
  #tax_tab.rick.unc.asvid,

  filt_tax_tab,
  chlor_tax_tab,
  MTunc_tax_tab,
  fileprefx,

  file = paste(fileprefx, "_allphyloseq.Rdata", sep = ""))


################################## no chloroplast here #######

### removed the chloroplast
fcount_tab = filt_fcount_tab


# replace 0 values with an estimate  ######### f.n0.count_tab
f.n0.count_tab <- cmultRepl(t(fcount_tab), method="CZM", label=0)

# generate the CLR values for plotting later
f.n0.clr <- codaSeq.clr(f.n0.count_tab)
####Returns a matrix of center log-ratio transformed data with samples by row.
##Each value is equivalent to log(x/gx) where gx is the geometic mean of the row vector X.
head(t(f.n0.clr), n=2)


###### reconstruct phyloseq again with the raw count and CZM method (replace 0 with an estimated).
#Transform into matrixes otu and tax tables
otu_mat <- as.matrix(fcount_tab)
tax_mat <- tax_tab[row.names(fcount_tab),]


OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(filt_sample_info_tab)  ### change to filter for ponds/tissues

