setwd("/home/mmorel/Documents/yabby/data/gymnosperms/")

Sys.setenv(
  PATH = paste(
    Sys.getenv("PATH"), "/home/mmorel/bin/miniconda3/bin", sep = ":"
  )
)

library(syntenet)
library(GenomicRanges)



aastringsetlist <- fasta2AAStringSetlist("gymno_synteny_data/")

annotation_cycas = read.csv("gymno_synteny_data/cycas.nosplicing.tsv", sep="\t") 
annotation_ginkgo = read.csv('gymno_synteny_data/Ginkgo_biloba.nosplicing.tsv', sep="\t") 
annotation_gnetum = read.csv('gymno_synteny_data/gnetum_montanum_nosplicing.tsv', sep="\t") 
annotation_sequoia = read.csv('gymno_synteny_data/sequoia.nosplicing.tsv', sep="\t") 
annotation_torreya = read.csv('gymno_synteny_data/torreya.nosplicing.tsv', sep="\t") 
annotation_welwitschia = read.csv('gymno_synteny_data/welwitschia.nosplicing.tsv', sep="\t") 

annot_cycas = makeGRangesFromDataFrame(annotation_cycas, seqnames.field="seqnames", start.field = "start", 
                                          end.field="end", strand.field="strand", keep.extra.columns=TRUE)
annot_ginkgo = makeGRangesFromDataFrame(annotation_ginkgo, seqnames.field = "seqnames", start.field = "start", 
                                     end.field = "end", strand.field = "strand", keep.extra.columns = TRUE)
annot_gnetum = makeGRangesFromDataFrame(annotation_gnetum, seqnames.field = "seqnames", start.field = "start", 
                                        end.field = "end", strand.field = "strand", keep.extra.columns = TRUE)
annot_sequoia = makeGRangesFromDataFrame(annotation_sequoia, seqnames.field = "seqnames", start.field = "start", 
                                        end.field = "end", strand.field = "strand", keep.extra.columns = TRUE)
annot_torreya = makeGRangesFromDataFrame(annotation_torreya, seqnames.field = "seqnames", start.field = "start", 
                                        end.field = "end", strand.field = "strand", keep.extra.columns = TRUE)
annot_welwitschia = makeGRangesFromDataFrame(annotation_welwitschia, seqnames.field = "seqnames", start.field = "start", 
                                        end.field = "end", strand.field = "strand", keep.extra.columns = TRUE)


annot_lst= list(annot_cycas, annot_ginkgo, annot_gnetum,annot_sequoia, annot_torreya, annot_welwitschia) #ordre alphabétique des diminutifs
names(annot_lst) <- c("cyc", "gink", "gnet","seq", "tor", "wel") #ordre alphabétique

names(aastringsetlist) <- c("cyc", "gink", "gnet", "seq", "tor", "wel") 
check_input(aastringsetlist, annot_lst)

#########

pdata <- process_input(aastringsetlist, annot_lst)

blast_list <- run_diamond(seq = pdata$seq, top_hits=5) #, evalue=0.001

write.table(do.call(rbind, blast_list), file='gymno_synteny_data/blast_table.tsv', quote=FALSE, sep='\t', row.names = FALSE)

######

'''
> names(blast_list)
"cyc_cyc"   "cyc_gink"  "cyc_gnet"  "gink_cyc"  "gink_gink" "gink_gnet" "gnet_cyc"  "gnet_gink" "gnet_gnet"

"cyc_cyc"   "cyc_gink"  "cyc_gnet"  "cyc_tor"   "cyc_wel"   "gink_cyc"  "gink_gink" "gink_gnet" "gink_tor"  "gink_wel"  "gnet_cyc"  "gnet_gink" "gnet_gnet" "gnet_tor" 
"gnet_wel"  "tor_cyc"   "tor_gink"  "tor_gnet"  "tor_tor"   "tor_wel"   "wel_cyc"   "wel_gink"  "wel_gnet"  "wel_tor"   "wel_wel" 

'''
list_A <- c()
a <- 0

for (i in 0:5) {
  for (j in 0:5) {
    a <- a + 1
    if (j > i) {
      list_A <- c(list_A, a)
    }
  }
}
print(list_A)

names(blast_list)


diamond_inter <- blast_list[list_A]
diamond_intra <- blast_list[c(1,8,15,22,29,36)]

value = 1e-5
anchors = 5
max_gaps = 25

dir.create(file.path('gymno_synteny_data/gymno_25_5_1e-5'))
intersyn <- interspecies_synteny(diamond_inter, pdata$annotation, inter_dir = "gymno_synteny_data/gymno_25_5_1e-5/" ,anchors=anchors, max_gaps=max_gaps, e_value=value)
#.collinearity
all_syn <- parse_collinearity(intersyn, as = "all")
write.table(all_syn , file='gymno_synteny_data/gymno_25_5_1e-5/inter_25_5_1e-5.tsv', quote=FALSE, sep='\t', row.names = FALSE)

####intrasynteny
dir.create(file.path('gymno_synteny_data/intra_gymno_25_5_1e-5'))
intrasyn <- intraspecies_synteny(diamond_intra, pdata$annotation,anchors=anchors, max_gaps=max_gaps, e_value=value, intra_dir = "gymno_synteny_data/intra_gymno_25_5_1e-5/" )
intra_syn <- parse_collinearity(intrasyn, as = "all")
write.table(intra_syn, file='gymno_synteny_data/intra_gymno_25_5_1e-5/intra_25_5_1e-5.tsv', quote=FALSE, sep='\t', row.names = FALSE)


value = 1e-3
anchors = 3
max_gaps = 55

dir.create(file.path('gymno_synteny_data/gymno_55_3_1e-3'))
intersyn <- interspecies_synteny(diamond_inter, pdata$annotation, inter_dir = "gymno_synteny_data/gymno_55_3_1e-3/" ,anchors=anchors, max_gaps=max_gaps, e_value=value)
#.collinearity

all_syn <- parse_collinearity(intersyn, as = "all")
write.table(all_syn , file='gymno_synteny_data/gymno_55_3_1e-3/inter_55_3_1e-3.tsv', quote=FALSE, sep='\t', row.names = FALSE)




####intrasynteny
dir.create(file.path('gymno_synteny_data/intra_gymno_55_5_1e-3'))
intrasyn <- intraspecies_synteny(diamond_intra, pdata$annotation,anchors=anchors, max_gaps=max_gaps, e_value=value, intra_dir = "gymno_synteny_data/intra_gymno_55_5_1e-3/" )
intra_syn <- parse_collinearity(intrasyn, as = "all")
write.table(intra_syn, file='gymno_synteny_data/intra_gymno_55_5_1e-3/intra_55_5_1e-3.tsv', quote=FALSE, sep='\t', row.names = FALSE)


