####Inter synteny and Intra synteny


##aristo vs ambo vs vitis
## intra synteny
##on the cluster
###conda activate syntenet_dev

library(syntenet)
library(GenomicRanges)

Sys.setenv(
    PATH = paste(
        Sys.getenv("PATH"), "/beegfs/data/mmorel/miniconda3/bin", sep = ":"
    )
)
setwd("/beegfs/data/mmorel/Infer_Anth/data/")

aastringsetlist <- fasta2AAStringSetlist("inter_species/")
annotation_ambo = read.csv("inter_species/amborella_trichopoda_nosplicing_63666.tsv", sep="\t")
annotation_aristo = read.csv("inter_species/aristolochia_fimbriata_nosplicing_61460.tsv", sep="\t")
annotation_vitis = read.csv("inter_species/vitis_vinifera_nosplicing_PN40024.tsv", sep="\t")


annot_ambo = makeGRangesFromDataFrame(annotation_ambo, seqnames.field="seqnames", start.field = "start", 
                         end.field="end", strand.field="strand", keep.extra.columns=TRUE)
annot_aristo = makeGRangesFromDataFrame(annotation_aristo, seqnames.field="seqnames", start.field = "start", 
                                 end.field="end", strand.field="strand", keep.extra.columns=TRUE)
annot_vitis = makeGRangesFromDataFrame(annotation_vitis, seqnames.field="seqnames", start.field = "start", 
                                 end.field="end", strand.field="strand", keep.extra.columns=TRUE)

annot_lst= list(annot_aristo, annot_ambo, annot_vitis) #ordre alphabétique des diminutifs

names(annot_lst) <- c("afi", "atr", "viv") #ordre alphabétique

names(aastringsetlist) <- c("afi", "atr", "viv")
check_input(aastringsetlist, annot_lst)


#########

pdata <- process_input(aastringsetlist, annot_lst)

blast_list <- run_diamond(seq = pdata$seq, top_hits=5)
write.table(do.call(rbind, blast_list), file='inter_species/ambo_aristo_vitis/blast_table.tsv', quote=FALSE, sep='\t', row.names = FALSE)

######


# > names(blast_list)
# [1] "afi_afi" "afi_atr" "afi_viv" "atr_afi" "atr_atr" "atr_viv" "viv_afi"
# [8] "viv_atr" "viv_viv"

diamond_intra <- blast_list[c(1,5,9)] #afi_afi, atr_atr, viv_viv

diamond_inter <- blast_list[c(2,3,6,8)] #afi_atr, afi_viv, atr_viv, viv_atr



value = 1e-3
anchors = 3
max_gaps = 55


intersyn <- interspecies_synteny(diamond_inter, pdata$annotation, inter_dir = "inter_species/ambo_aristo_vitis/55_3_1e-3/" ,anchors=anchors, max_gaps=max_gaps, e_value=value)

intrasyn <- intraspecies_synteny(diamond_intra, pdata$annotation, intra_dir = "inter_species/ambo_aristo_vitis/55_3_1e-3/", anchors=anchors, max_gaps=max_gaps, e_value=value)


#.collinearity

all_syn <- parse_collinearity(intersyn, as = "all")
write.table(all_syn , file='inter_species/ambo_aristo_nymphaea_vitis/55_3_1e-3/inter_55_3_1e-3.tsv', quote=FALSE, sep='\t', row.names = FALSE)





value = 5e-2
anchors = 3
max_gaps = 55


intersyn <- interspecies_synteny(diamond_inter, pdata$annotation, inter_dir = "inter_species/ambo_aristo_vitis/55_3_5e-2/" ,anchors=anchors, max_gaps=max_gaps, e_value=value)
intrasyn <- intraspecies_synteny(diamond_intra, pdata$annotation, intra_dir = "inter_species/ambo_aristo_vitis/55_3_5e-2/", anchors=anchors, max_gaps=max_gaps, e_value=value)


#.collinearity

all_syn <- parse_collinearity(intersyn, as = "all")
write.table(all_syn , file='inter_species/ambo_aristo_nymphaea_vitis/55_3_1e-3/inter_55_3_1e-3.tsv', quote=FALSE, sep='\t', row.names = FALSE)





#yabby_afi = list("afi_afi_Af03G134500", "afi_afi_Af01G026600", "afi_afi_Af02G072800", "afi_afi_Af05G185300","afi_afi_Af07G024800")
#INO, YAB5, FIL, YAB2, CRC
#yabby_atr = list("atr_atr_94506842", "atr_atr_94491478", "atr_atr_94503013", "atr_atr_94473924", "atr_atr_94507018")
#YAB2, CRC, INO, YAB5, FIL
#yabby_viv = list("viv_Vitvi06g00972", "viv_Vitvi02g00510", "viv_Vitvi15g00708", "viv_Vitvi11g00492", "viv_Vitvi01g00013",  "viv_Vitvi08g00274", "viv_viv_Vitvi01g00703")
# YAB2 FIL FIL YAB5 CRC YAB2 INO
#yabby_cycas = list("cyp_CYCAS_007102", "cyp_CYCAS_008411", "cyp_CYCAS_022386")