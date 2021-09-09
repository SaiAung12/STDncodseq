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
library(ggsci)
library(scales)
library(tidyr)
library(dplyr)
library(compositions)

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

### READING IN OUR DATA
setwd("/home/sylvester/test/PG01-D4-STDncodseq04-01")

##output directory
outdir = "/home/sylvester/test/PG01-D4-STDncodseq04-01/dada2fspc1/"

if (file.exists(outdir)){
    setwd(file.path(outdir))
} else {
    dir.create(file.path(outdir))
    setwd(file.path(outdir))
}


load("/home/pvannameihome/PG01-D4-16S-STDncodseq03-01/dada2f230/202102-pg01-D4-16S-rmchloMTunc-Q1f230-STDnCoda01-minprop-001minOc02OTUs1110_aldextest.Rdata")
load("/home/pvannameihome/PG01-D4-16S-STDncodseq03-01/dada2f230/202102-pg01-D4-16S-rmchloMTunc-Q1f230-STDnCoda01-minprop-001minOc02OTUs1110_allphyloseq.Rdata")

testcolk = pal_simpsons("springfield", alpha = 0.75)(16)

wkcolpal = testcolk[c(1, 5, 6,
                      15, 2, 3, 9, 8,
                      11, 10, 16, 13, 12,
                      7, 4, 14)]

names(wkcolpal) = c("n0", "n4", "n6", 0:12)

### INPUT should be from
###  STDncodseq03  ###
##### set params
fpond = "pg01-D4"
ftiss = "16S-rmchloMTunc"
prop = "Q1f230-STDnCoda01-minprop"
prop1 = "Q1f230-STDnCoda01-cutoffs"

fminprop = "001"
fminoc =  "02"
################################## no chloroplast and MTunc here #######
### removed the chloroplast
#remove chloroplast  so change the fcount_tab to
#reset the number of OTUs after removed ### 752 OTUs
nbrOTUs = dim(fcount_tab)[1]
nbrOTUs

#1110

fileprefx = paste("202102-", fpond, "-",
                  ftiss, "-",
                  "Q1f275-STDnCoda01-minprop", fminprop, "minOc", fminoc, "OTUs", nbrOTUs, sep = "")
fileprefx


## print tax
tax_table(f2relCZM.PG01p3Q1_noChloMTunc_physeq)
## f2relCZM.PG01p3Q1_noChloMTunc_physeq  is the non - Chloroplast and MTunc


###############

####### functions  ###################################

find.top.asvID <- function(x,num){
  require(phyloseq)
  require(magrittr)
  totu <- otu_table(x)
  j1 <- apply(totu,2,sort,index.return=T, decreasing=T) # modifying which.max to return a list of sorted index
  j2 <- lapply(j1,'[[',"ix") # select for index
  topn = j2 %>%  lapply(head,n=num)
  mid = NA
  if(length(unique(unlist(topn))) == 1 ){
    mid = row.names(totu@.Data[c(unique(unlist(topn)), unique(unlist(topn))), 1:2])
    mid = mid[1]
  } else {
    mid = row.names(totu@.Data[unique(unlist(topn)), ])
  }
  return(mid)

}

## get values
get.top.asvValue <- function(x,num){
  require(phyloseq)
  require(magrittr)
  totu <- otu_table(x)
  motu = as(totu, "matrix")
  j1 <- apply(totu,2,sort,index.return=T, decreasing=T) # modifying which.max to return a list of sorted index
  j2 <- lapply(j1,'[[',"ix") # select for index
  topn = j2 %>%  lapply(head,n=num)

  vtop = lapply(names(topn), function(n)  { motu[topn[[n]], n ] })

  vdf <- data.frame(matrix(unlist(vtop), ncol=length(vtop), byrow=F))
  names(vdf) <- names(topn)
  return(vdf)

}


### from PG01DS0102

## change string NA to NA
bktax01 = tax_table(f2relCZM.PG01p3Q1_noChloMTunc_physeq)
tmptaxtab = tax_table(f2relCZM.PG01p3Q1_noChloMTunc_physeq)
head(tmptaxtab@.Data[,7])
for (ip in 1:7) {
  tmptaxtab@.Data[which(tmptaxtab@.Data[,ip] == "NA"), ip] = NA
}
tax_table(f2relCZM.PG01p3Q1_noChloMTunc_physeq) = tmptaxtab

####

### get best hit at with microbiomeutilities::format_to_besthit
f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit = microbiomeutilities::format_to_besthit(f2relCZM.PG01p3Q1_noChloMTunc_physeq)
head(tax_table(f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit))
head(otu_table(f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit))
## count at species level
# number of ASV

taxtable.besthit = tax_table(f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit)

bestspc.cnt = table(taxtable.besthit@.Data[,7])
names(bestspc.cnt)
tax.kcol = which(grepl("k__",names(bestspc.cnt)))
tax.pcol = which(grepl("p__",names(bestspc.cnt)))
tax.ccol = which(grepl("c__",names(bestspc.cnt)))
tax.ocol = which(grepl("o__",names(bestspc.cnt)))
tax.fcol = which(grepl("f__",names(bestspc.cnt)))
tax.gcol = which(grepl("g__",names(bestspc.cnt)))
tax.scol = which(!grepl("__",names(bestspc.cnt)))

bestspc.tab = c(sum(bestspc.cnt[tax.kcol]),
                sum(bestspc.cnt[tax.pcol]),
                sum(bestspc.cnt[tax.ccol]),
                sum(bestspc.cnt[tax.ocol]),
                sum(bestspc.cnt[tax.fcol]),
                sum(bestspc.cnt[tax.gcol]),
                sum(bestspc.cnt[tax.scol])
)
names(bestspc.tab) = colnames(taxtable.besthit)[-8]

write.table(bestspc.tab, file =paste(fileprefx, "-cntnbr-NoChloMTunc_ASV_atbestTaxRank.tab", sep = ""), sep = "\t" )


################################ at ASV levels ####
f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf = phy_to_ldf(f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit, transform.counts = NULL)
f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf$nID = f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf$OTUID

f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf$wk = factor(f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf$wk)
## set new ordering by ponds and by week
f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf$pdwk = paste(f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf$pd,
                                                         f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf$wk,
                                                         sep = "-")


for (iq in which(!grepl("__", f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf$Species))) {
  f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf$nID[iq] = paste(sep = ":", (strsplit(
    as.character(f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf$OTUID[iq]),
    ":", fixed = TRUE))[[1]],
    f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf$Species[iq])

}

f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf$oID = f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf$OTUID
for (iq in 1:dim(f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf)) {
  f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf$oID[iq] = (strsplit(
    as.character(f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf$OTUID[iq]),
    ":", fixed = TRUE))[[1]]
  f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf$oID[iq] = gsub("^OTU-", "",
                                                             f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf$oID[iq])
}


## CHANGE NAME FOR ordering
#taxtable.besthit.mtrx = as(taxtable.besthit, "matrix")
head(f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf)
taxtable.besthit.mtrx.lst = as.vector(f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf$best_hit)
taxtable.besthit.spp.lst = as.vector(f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf$Species)
names(taxtable.besthit.mtrx.lst) = NULL
names(taxtable.besthit.spp.lst) = NULL

tax.kcol = which(grepl("k__",taxtable.besthit.mtrx.lst))
tax.pcol = which(grepl("p__",taxtable.besthit.mtrx.lst))
tax.ccol = which(grepl("c__",taxtable.besthit.mtrx.lst))
tax.ocol = which(grepl("o__",taxtable.besthit.mtrx.lst))
tax.fcol = which(grepl("f__",taxtable.besthit.mtrx.lst))
tax.gcol = which(grepl("g__",taxtable.besthit.spp.lst))
tax.scol = which(!grepl("__",taxtable.besthit.spp.lst))

str(taxtable.besthit.mtrx.lst)
namesplit = strsplit(taxtable.besthit.mtrx.lst, ":", fixed = TRUE)
str(namesplit)

hn = namesplit %>%  lapply(head,n=1)
tn = namesplit %>%  lapply(tail,n=1)
tn = unlist(tn)
hn = unlist(hn)
tnn = strsplit(tn, "__", fixed = TRUE)
tnf = unlist(tnn %>%  lapply(tail,n=1))

nnres = data.frame(f=taxtable.besthit.mtrx.lst, hn=hn, tn=tn, tnf=tnf, sp = taxtable.besthit.spp.lst, nID=f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf$nID,
                   check.names = F)
head(nnres)
head(nnres[tax.scol, ] )
dim(nnres)
slx.sp = nnres[tax.scol, ]
slx.sp.v = as.vector(slx.sp[order(slx.sp$sp), 6 ])
slx.sp.uniqID = unique(slx.sp.v, fromLast = FALSE)
head(slx.sp.uniqID)

####
slx.sp = nnres[tax.kcol, ]
slx.sp.v = as.vector(slx.sp[order(slx.sp$tnf), 1 ])
slx.k.uniqID = unique(slx.sp.v, fromLast = FALSE)
head(slx.k.uniqID)
slx.sp = nnres[tax.pcol, ]
slx.sp.v = as.vector(slx.sp[order(slx.sp$tnf), 1 ])
slx.p.uniqID = unique(slx.sp.v, fromLast = FALSE)
head(slx.p.uniqID)
slx.sp = nnres[tax.ccol, ]
slx.sp.v = as.vector(slx.sp[order(slx.sp$tnf), 1 ])
slx.c.uniqID = unique(slx.sp.v, fromLast = FALSE)
head(slx.c.uniqID)
slx.sp = nnres[tax.ocol, ]
slx.sp.v = as.vector(slx.sp[order(slx.sp$tnf), 1 ])
slx.o.uniqID = unique(slx.sp.v, fromLast = FALSE)
head(slx.o.uniqID)
slx.sp = nnres[tax.fcol, ]
slx.sp.v = as.vector(slx.sp[order(slx.sp$tnf), 1 ])
slx.f.uniqID = unique(slx.sp.v, fromLast = FALSE)
head(slx.f.uniqID)
slx.sp = nnres[tax.gcol, ]
slx.sp.v = as.vector(slx.sp[order(slx.sp$tnf), 1 ])
slx.g.uniqID = unique(slx.sp.v, fromLast = FALSE)
head(slx.g.uniqID)


#sort(nnres[tax.scol, ] , index.return=T)


f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf$nID <- factor(f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf$nID,
                                                          levels = c(slx.sp.uniqID,
                                                                     slx.g.uniqID,
                                                                     slx.f.uniqID,
                                                                     slx.o.uniqID,
                                                                     slx.c.uniqID,
                                                                     slx.p.uniqID,
                                                                     slx.k.uniqID

                                                          ))

head(f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf)
f2rel.n2oID = unique(select(f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf, nID, oID,
                            Domain, Phylum, Class, Order, Family, Genus,  Species))


## forgraphic below
x.sampledat = sample_data(f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit)
x.sampledat$pdwk=paste(x.sampledat$pd, x.sampledat$wk, sep = "-")
pdwk.order = unique(x.sampledat$pdwk)
samp.by.pdwk = row.names(x.sampledat[order(x.sampledat$pdwk), ])



### save objects for besthit taxonomic levels at ASVs
save(
  f2relCZM.PG01p3Q1_noChloMTunc_physeq,
  f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit,
  f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf,

  file = paste(fileprefx, "_besthits_seq04_phyloseq.Rdata", sep = ""))

iq =3

grpName= paste("top", "-", as.character(iq),"_ASV", sep = "")

top_taxa_id = find.top.asvID(f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit, iq)
#### top taxa
#top_taxa_id = find.top.asvID(f2relCZM.wfs0102physeq.besthit, 1)
#find.top.asvID(wfs0102physeq, 2)
#find.top.asvID(wfs0102physeq, 5)
# filter(psmelt(f2relCZM.Nchlo_pg01ds0102Q1physeq.besthit),
#        OTU %in% top_taxa_id)


p.box <- ggstripchart(filter(f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf, OTUID %in% top_taxa_id),
                      "Sam_rep", "Abundance",
                      shape = "source",
                      facet.by = "nID", color = "wk",
                      palette = wkcolpal,
                      title = grpName,
                      order = samp.by.pdwk
)

datasub = filter(f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf, OTUID %in% top_taxa_id)
datasub$nameg =datasub$OTUID

for (iz in 1:(dim(datasub)[1])){
  cnID = gsub("^OTU-", "", datasub$nID[iz])
  datasub$nameg[iz] = paste(collapse = "\n", ((strsplit(as.character(cnID), ":", fixed = TRUE))[[1]]))
}

nlabel = datasub$nameg
names(nlabel) = datasub$nID


p.box + rremove("x.text") + scale_y_continuous(breaks=c(0,0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 1.0), trans = "log10") +
  facet_wrap(~nID, ncol = 3,
             labeller = labeller(nID=nlabel)) +
  #  facet_wrap(~nID, ncol = 3) +
  theme(
    strip.text.x = element_text(size = 8),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "grey"),
    panel.grid.major.x = element_blank()
  )


ggsave(filename=paste(fileprefx, "-", grpName, "NoChloro_MTunc_ASV_abundance_sc.pdf", sep = ""), width=210, height=1500, units="mm",
       limitsize = FALSE)


topx.ldf = filter(f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf, OTUID %in% top_taxa_id)
names(topx.ldf)
topx.ldf.s1 = topx.ldf[c(1, 21, 10, 11)]
topx_wide <- spread(topx.ldf.s1, Sam_rep, Abundance)
topx_wide[1:2, ]

row.names(topx_wide) = topx_wide$oID
topx.wide.s1 = merge(topx_wide, (data.frame(as(otu_table(f2.PG01p3Q1_noChloMTunc_physeq), "matrix")))[topx_wide$oID, ],
                     by =0)

row.names(topx.wide.s1)  = topx.wide.s1$Row.names
topx.wide.s1 = topx.wide.s1[-1]
topx.wide.s1.tax = merge(topx.wide.s1,
                         data.frame(tax_table(f2.PG01p3Q1_noChloMTunc_physeq)[topx_wide$oID, ]),
                         by = 0)
row.names(topx.wide.s1.tax) = topx.wide.s1.tax$Row.names
topx.wide.s1.tax = topx.wide.s1.tax[-1]
topx.wide.s1.tax[1:2, ]
write.table(topx.wide.s1.tax, file =paste(fileprefx, "-", grpName, "-NoChloro_MTunc_ASV.tab", sep = ""), sep = "\t" )


top3_taxa_id = top_taxa_id

iq =5

grpName= paste("top", "-", as.character(iq),"_ASV", sep = "")

top5_taxa_id = find.top.asvID(f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit, iq)
top_taxa_id = subset(top5_taxa_id, !(top5_taxa_id %in% top3_taxa_id))


p.box <- ggstripchart(filter(f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf, OTUID %in% top_taxa_id),
                      "Sam_rep", "Abundance",
                      shape = "source",
                      facet.by = "nID", color = "wk",
                      palette = wkcolpal,
                      title = grpName,
                      order = samp.by.pdwk

)

datasub = filter(f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf, OTUID %in% top_taxa_id)
datasub$nameg =datasub$OTUID

for (iz in 1:(dim(datasub)[1])){
  cnID = gsub("^OTU-", "", datasub$nID[iz])
  datasub$nameg[iz] = paste(collapse = "\n", ((strsplit(as.character(cnID), ":", fixed = TRUE))[[1]]))
}

nlabel = datasub$nameg
names(nlabel) = datasub$nID


p.box + rremove("x.text") + scale_y_continuous(breaks=c(0,0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 1.0), trans = "log10") +
  facet_wrap(~nID, ncol = 3,
             labeller = labeller(nID=nlabel)) +
  #  facet_wrap(~nID, ncol = 3) +
  theme(
    strip.text.x = element_text(size = 8),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "grey"),
    panel.grid.major.x = element_blank()
  )


ggsave(filename=paste(fileprefx, "-", grpName, "NoChloro_MTunc_ASV_abundance_sc.pdf", sep = ""), width=210, height=900, units="mm",
       limitsize = FALSE)


topx.ldf = filter(f2relCZM.PG01p3Q1_noChloMTunc_physeq.besthit.ldf, OTUID %in% top_taxa_id)
names(topx.ldf)
topx.ldf.s1 = topx.ldf[c(1, 22, 10, 11)]
topx_wide <- spread(topx.ldf.s1, Sam_rep, Abundance)
topx_wide[1:2, ]

row.names(topx_wide) = topx_wide$oID
topx.wide.s1 = merge(topx_wide, (data.frame(as(otu_table(f2.PG01p3Q1_noChloMTunc_physeq), "matrix")))[topx_wide$oID, ],
                     by =0)

row.names(topx.wide.s1)  = topx.wide.s1$Row.names
topx.wide.s1 = topx.wide.s1[-1]
topx.wide.s1.tax = merge(topx.wide.s1,
                         data.frame(tax_table(f2.PG01p3Q1_noChloMTunc_physeq)[topx_wide$oID, ]),
                         by = 0)
row.names(topx.wide.s1.tax) = topx.wide.s1.tax$Row.names
topx.wide.s1.tax = topx.wide.s1.tax[-1]
topx.wide.s1.tax[1:2, ]
write.table(topx.wide.s1.tax, file =paste(fileprefx, "-", grpName, "-NoChloro_MTunc_ASV.tab", sep = ""), sep = "\t" )

