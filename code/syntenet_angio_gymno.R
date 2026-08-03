setwd("/home/mmorel/Documents/yabby/data/")

Sys.setenv(
  PATH = paste(
    Sys.getenv("PATH"), "/home/mmorel/bin/miniconda3/bin", sep = ":"
  )
)

library(syntenet)
library(GenomicRanges)


aastringsetlist <- fasta2AAStringSetlist("synteny_angio_gymno/")

annotation_ambo = read.csv("synteny_angio_gymno/amborella_trichopoda_nosplicing_63666.tsv", sep="\t")
annotation_cycas = read.csv("synteny_angio_gymno/cycas.nosplicing.tsv", sep="\t")
annotation_ginkgo = read.csv('synteny_angio_gymno/Ginkgo_biloba.nosplicing.tsv', sep="\t") 
annotation_illicium = read.csv("synteny_angio_gymno/illicium_verum_no_splicing.tsv", sep="\t")
annotation_nymphaea = read.csv("synteny_angio_gymno/nymphaea_colorata.tsv", sep="\t")
#annotation_sequoia = read.csv('synteny_angio_gymno/sequoia.nosplicing.tsv', sep="\t") 
#annotation_vitis = read.csv("synteny_angio_gymno/vitis_vinifera_nosplicing_PN40024.tsv", sep="\t")


annot_ambo = makeGRangesFromDataFrame(annotation_ambo, seqnames.field="seqnames", start.field = "start", 
                                      end.field="end", strand.field="strand", keep.extra.columns=TRUE)

annot_cycas = makeGRangesFromDataFrame(annotation_cycas, seqnames.field="seqnames", start.field = "start", 
                                      end.field="end", strand.field="strand", keep.extra.columns=TRUE)

annot_ginkgo = makeGRangesFromDataFrame(annotation_ginkgo, seqnames.field = "seqnames", start.field = "start", 
                                        end.field = "end", strand.field = "strand", keep.extra.columns = TRUE)
annot_illicium = makeGRangesFromDataFrame(annotation_illicium, seqnames.field = "seqnames", start.field = "start", 
                                          end.field = "end", strand.field = "strand", keep.extra.columns = TRUE)
annot_nymphaea = makeGRangesFromDataFrame(annotation_nymphaea, seqnames.field="seqnames", start.field = "start", 
                                          end.field="end", strand.field="strand", keep.extra.columns=TRUE)
#annot_sequoia = makeGRangesFromDataFrame(annotation_sequoia, seqnames.field = "seqnames", start.field = "start", 
#                                         end.field = "end", strand.field = "strand", keep.extra.columns = TRUE)
#annot_vitis = makeGRangesFromDataFrame(annotation_vitis, seqnames.field="seqnames", start.field = "start", 
#                                       end.field="end", strand.field="strand", keep.extra.columns=TRUE)




#annot_lst= list(annot_ambo, annot_ginkgo, annot_illicium, annot_nymphaea, annot_sequoia, annot_vitis) #ordre alphabétique des diminutifs
annot_lst= list(annot_ambo, annot_cycas, annot_ginkgo, annot_illicium, annot_nymphaea) #ordre alphabétique des diminutifs

#names(annot_lst) <- c("atr", "gink","ive", "nco", "seq", "viv") #ordre alphabétique
names(annot_lst) <- c("atr", "cyc", "gink","ive", "nco") #ordre alphabétique

#names(aastringsetlist) <- c("atr", "gink","ive", "nco", "seq", "viv") 
names(aastringsetlist) <-c("atr", "cyc", "gink","ive", "nco")
check_input(aastringsetlist, annot_lst)

#########

pdata <- process_input(aastringsetlist, annot_lst)

blast_list <- run_diamond(seq = pdata$seq, top_hits=5) #, evalue=0.001

write.table(do.call(rbind, blast_list), file='synteny_angio_gymno/atr_cyc_gink_ive_nco/blast_table.tsv', quote=FALSE, sep='\t', row.names = FALSE)

######

#> names(blast_list)
#"cyc_cyc"   "cyc_gink"  "cyc_gnet"  "gink_cyc"  "gink_gink" "gink_gnet" "gnet_cyc"  "gnet_gink" "gnet_gnet"

#"cyc_cyc"   "cyc_gink"  "cyc_gnet"  "cyc_tor"   "cyc_wel"   "gink_cyc"  "gink_gink" "gink_gnet" "gink_tor"  "gink_wel"  "gnet_cyc"  "gnet_gink" "gnet_gnet" "gnet_tor" 
#"gnet_wel"  "tor_cyc"   "tor_gink"  "tor_gnet"  "tor_tor"   "tor_wel"   "wel_cyc"   "wel_gink"  "wel_gnet"  "wel_tor"   "wel_wel" 

list_A <- c()
a <- 0

for (i in 0:4) {
  for (j in 0:4) {
    a <- a + 1
    if (j > i) {
      list_A <- c(list_A, a)
    }
  }
}
print(list_A)

names(blast_list)


diamond_inter <- blast_list[list_A]
diamond_intra <- blast_list[c(1,7,13,19,25)]

value = 1e-5
anchors = 5
max_gaps = 25

dir.create(file.path('synteny_angio_gymno/angio_gymno_25_5_1e-5'))
intersyn <- interspecies_synteny(diamond_inter, pdata$annotation, inter_dir = "synteny_angio_gymno/angio_gymno_25_5_1e-5/",
                                 anchors=anchors, max_gaps=max_gaps, e_value=value)
#.collinearity
all_syn <- parse_collinearity(intersyn, as = "all")
write.table(all_syn , file='synteny_angio_gymno/angio_gymno_25_5_1e-5/inter_25_5_1e-5.tsv', quote=FALSE, sep='\t', row.names = FALSE)

####intrasynteny
dir.create(file.path('synteny_angio_gymno/intra_angio_gymno_25_5_1e-5'))
intrasyn <- intraspecies_synteny(diamond_intra, pdata$annotation,anchors=anchors, max_gaps=max_gaps, e_value=value,
                                 intra_dir = "synteny_angio_gymno/intra_angio_gymno_25_5_1e-5/" )
intra_syn <- parse_collinearity(intrasyn, as = "all")
write.table(intra_syn, file='synteny_angio_gymno/intra_angio_gymno_25_5_1e-5/intra_25_5_1e-5.tsv', quote=FALSE, sep='\t', row.names = FALSE)


value = 1e-3
anchors = 3
max_gaps = 55

dir.create(file.path('synteny_angio_gymno/angio_gymno_55_3_1e-3'))
intersyn <- interspecies_synteny(diamond_inter, pdata$annotation, inter_dir = "synteny_angio_gymno/angio_gymno_55_3_1e-3/",
                                 anchors=anchors, max_gaps=max_gaps, e_value=value)
#.collinearity

all_syn <- parse_collinearity(intersyn, as = "all")
write.table(all_syn , file='synteny_angio_gymno/angio_gymno_55_3_1e-3/inter_55_3_1e-3.tsv', quote=FALSE, sep='\t', row.names = FALSE)


####intrasynteny
dir.create(file.path('synteny_angio_gymno/intra_angio_gymno_55_5_1e-3'))
intrasyn <- intraspecies_synteny(diamond_intra, pdata$annotation,anchors=anchors, max_gaps=max_gaps,
                                 e_value=value, intra_dir = "synteny_angio_gymno/intra_angio_gymno_55_5_1e-3/" )
intra_syn <- parse_collinearity(intrasyn, as = "all")
write.table(intra_syn, file='synteny_angio_gymno/intra_angio_gymno_55_5_1e-3/intra_55_5_1e-3.tsv', quote=FALSE, sep='\t', row.names = FALSE)
