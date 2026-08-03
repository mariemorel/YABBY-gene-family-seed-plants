########################### intraspecies illicium verum

###conda activate syntenet_dev
setwd("/beegfs/data/mmorel/Infer_Anth/data/illicium/")

Sys.setenv(
    PATH = paste(
        Sys.getenv("PATH"), "/beegfs/data/mmorel/miniconda3/bin", sep = ":"
    )
)

library(syntenet)
library(GenomicRanges)


aastringsetlist <- fasta2AAStringSetlist("data_ive/")
annotation_ive = read.csv("data_ive/illicium_verum_no_splicing.tsv", sep="\t")
annot_ive= makeGRangesFromDataFrame(annotation_ive, seqnames.field="seqnames", start.field = "start", 
                                 end.field="end", strand.field="strand", keep.extra.columns=TRUE)

annot_lst= list(annot_ive)
names(annot_lst) <- c("ive")
names(aastringsetlist) <- c("ive")

check_input(aastringsetlist, annot_lst)

###################

pdata <- process_input(aastringsetlist, annot_lst)
blast_list <- run_diamond(seq = pdata$seq, top_hits=5)

write.table(blast_list, file='data_ive/blast_table.tsv', quote=FALSE, sep='\t', row.names = FALSE)

#####################

value = 1e-5
anchors = 5
max_gaps = 25

intrasyn <- intraspecies_synteny(blast_list, pdata$annotation,anchors=anchors, max_gaps=max_gaps, e_value=value, intra_dir = "data_ive/ive_25_5_1e-5/" )

ive_syn <- parse_collinearity(intrasyn, as = "all")
write.table(ive_syn, file='data_ive/ive_25_5_1e-5/intra_25_5_1e-5.tsv', quote=FALSE, sep='\t', row.names = FALSE)


value = 1e-3
anchors = 3
max_gaps = 55

intrasyn <- intraspecies_synteny(blast_list, pdata$annotation,anchors=anchors, max_gaps=max_gaps, e_value=value, intra_dir = "data_ive/ive_55_3_1e-3/" )

ive_syn <- parse_collinearity(intrasyn, as = "all")
write.table(ive_syn, file='data_ive/ive_55_3_1e-3/intra_55_3_1e-3.tsv', quote=FALSE, sep='\t', row.names = FALSE)

##########intraspecies syntenet_dev aristolochia contorta

###conda activate syntenet_dev
setwd("/beegfs/data/mmorel/Infer_Anth/data/aristo_contorta/")

Sys.setenv(
    PATH = paste(
        Sys.getenv("PATH"), "/beegfs/data/mmorel/miniconda3/bin", sep = ":"
    )
)

library(syntenet)
library(GenomicRanges)


aastringsetlist <- fasta2AAStringSetlist("data_aco/")
annotation_aco = read.csv("data_aco/aristo_contorta_no_splicing.tsv", sep="\t")
annot_aco= makeGRangesFromDataFrame(annotation_aco, seqnames.field="seqnames", start.field = "start", 
                                 end.field="end", strand.field="strand", keep.extra.columns=TRUE)

annot_lst= list(annot_aco)
names(annot_lst) <- c("aco")
names(aastringsetlist) <- c("aco")

check_input(aastringsetlist, annot_lst)

###################

pdata <- process_input(aastringsetlist, annot_lst)
blast_list <- run_diamond(seq = pdata$seq, top_hits=5)

write.table(blast_list, file='data_aco/blast_table.tsv', quote=FALSE, sep='\t', row.names = FALSE)

#####################

value = 1e-5
anchors = 5
max_gaps = 25

intrasyn <- intraspecies_synteny(blast_list, pdata$annotation,anchors=anchors, max_gaps=max_gaps, e_value=value, intra_dir = "data_aco/aco_25_5_1e-5/" )

aco_syn <- parse_collinearity(intrasyn, as = "all")
write.table(aco_syn, file='data_aco/aco_25_5_1e-5/intra_25_5_1e-5.tsv', quote=FALSE, sep='\t', row.names = FALSE)



value = 1e-5
anchors = 3
max_gaps = 35

intrasyn <- intraspecies_synteny(blast_list, pdata$annotation,anchors=anchors, max_gaps=max_gaps, e_value=value, intra_dir = "data_aco/aco_35_3_1e-5/" )

aco_syn <- parse_collinearity(intrasyn, as = "all")
write.table(aco_syn, file='data_aco/aco_35_3_1e-5/intra_35_3_1e-5.tsv', quote=FALSE, sep='\t', row.names = FALSE)


value = 1e-5
anchors = 3
max_gaps = 50

intrasyn <- intraspecies_synteny(blast_list, pdata$annotation,anchors=anchors, max_gaps=max_gaps, e_value=value, intra_dir = "data_aco/aco_50_3_1e-5/" )

aco_syn <- parse_collinearity(intrasyn, as = "all")
write.table(aco_syn, file='data_aco/aco_50_3_1e-5/intra_50_3_1e-5.tsv', quote=FALSE, sep='\t', row.names = FALSE)


value = 1e-3
anchors = 3
max_gaps = 35

intrasyn <- intraspecies_synteny(blast_list, pdata$annotation,anchors=anchors, max_gaps=max_gaps, e_value=value, intra_dir = "data_aco/aco_35_3_1e-3/" )

aco_syn <- parse_collinearity(intrasyn, as = "all")
write.table(aco_syn, file='data_aco/aco_35_3_1e-3/intra_35_3_1e-3.tsv', quote=FALSE, sep='\t', row.names = FALSE)


value = 1e-3
anchors = 3
max_gaps = 55

intrasyn <- intraspecies_synteny(blast_list, pdata$annotation,anchors=anchors, max_gaps=max_gaps, e_value=value, intra_dir = "data_aco/aco_55_3_1e-3/" )

aco_syn <- parse_collinearity(intrasyn, as = "all")
write.table(aco_syn, file='data_aco/aco_55_3_1e-3/intra_55_3_1e-3.tsv', quote=FALSE, sep='\t', row.names = FALSE)






######################annona


setwd("/beegfs/data/mmorel/Infer_Anth/data/annona_cherimola/")


aastringsetlist <- fasta2AAStringSetlist("data_ach/")
annotation_ach = read.csv("data_ach/annona_cherimola_no_splicing.tsv", sep="\t")
annot_ach= makeGRangesFromDataFrame(annotation_ach, seqnames.field="seqnames", start.field = "start", 
                                 end.field="end", strand.field="strand", keep.extra.columns=TRUE)

annot_lst= list(annot_ach)
names(annot_lst) <- c("ach")
names(aastringsetlist) <- c("ach")

check_input(aastringsetlist, annot_lst)

###################

pdata <- process_input(aastringsetlist, annot_lst)
blast_list <- run_diamond(seq = pdata$seq, top_hits=5)

write.table(blast_list, file='data_ach/blast_table.tsv', quote=FALSE, sep='\t', row.names = FALSE)

#####################

value = 1e-5
anchors = 5
max_gaps = 25

intrasyn <- intraspecies_synteny(blast_list, pdata$annotation,anchors=anchors, max_gaps=max_gaps, e_value=value, intra_dir = "data_ach/ach_25_5_1e-5/" )

ach_syn <- parse_collinearity(intrasyn, as = "all")
write.table(ach_syn, file='data_ach/ach_25_5_1e-5/intra_25_5_1e-5.tsv', quote=FALSE, sep='\t', row.names = FALSE)



value = 1e-5
anchors = 3
max_gaps = 35

intrasyn <- intraspecies_synteny(blast_list, pdata$annotation,anchors=anchors, max_gaps=max_gaps, e_value=value, intra_dir = "data_ach/ach_35_3_1e-5/" )

ach_syn <- parse_collinearity(intrasyn, as = "all")
write.table(ach_syn, file='data_ach/ach_35_3_1e-5/intra_35_3_1e-5.tsv', quote=FALSE, sep='\t', row.names = FALSE)


value = 1e-5
anchors = 3
max_gaps = 50

intrasyn <- intraspecies_synteny(blast_list, pdata$annotation,anchors=anchors, max_gaps=max_gaps, e_value=value, intra_dir = "data_ach/ach_50_3_1e-5/" )

ach_syn <- parse_collinearity(intrasyn, as = "all")
write.table(ach_syn, file='data_ach/ach_50_3_1e-5/intra_50_3_1e-5.tsv', quote=FALSE, sep='\t', row.names = FALSE)


value = 1e-3
anchors = 3
max_gaps = 35

intrasyn <- intraspecies_synteny(blast_list, pdata$annotation,anchors=anchors, max_gaps=max_gaps, e_value=value, intra_dir = "data_ach/ach_35_3_1e-3/" )

ach_syn <- parse_collinearity(intrasyn, as = "all")
write.table(ach_syn, file='data_ach/ach_35_3_1e-3/intra_35_3_1e-3.tsv', quote=FALSE, sep='\t', row.names = FALSE)


value = 1e-3
anchors = 3
max_gaps = 55

intrasyn <- intraspecies_synteny(blast_list, pdata$annotation,anchors=anchors, max_gaps=max_gaps, e_value=value, intra_dir = "data_ach/ach_55_3_1e-3/" )

ach_syn <- parse_collinearity(intrasyn, as = "all")
write.table(ach_syn, file='data_ach/ach_55_3_1e-3/intra_55_3_1e-3.tsv', quote=FALSE, sep='\t', row.names = FALSE)



#################nymphaea colorata


setwd("/beegfs/data/mmorel/Infer_Anth/data/colorata/")


aastringsetlist <- fasta2AAStringSetlist("data_nco/")
annotation_nco = read.csv("data_nco/nymphaea_colorata.tsv", sep="\t")
annot_nco= makeGRangesFromDataFrame(annotation_nco, seqnames.field="seqnames", start.field = "start", 
                                 end.field="end", strand.field="strand", keep.extra.columns=TRUE)

annot_lst= list(annot_nco)
names(annot_lst) <- c("nco")
names(aastringsetlist) <- c("nco")

check_input(aastringsetlist, annot_lst)

###################

pdata <- process_input(aastringsetlist, annot_lst)
blast_list <- run_diamond(seq = pdata$seq, top_hits=5)

write.table(blast_list, file='data_nco/blast_table.tsv', quote=FALSE, sep='\t', row.names = FALSE)

#####################

value = 1e-5
anchors = 5
max_gaps = 25

intrasyn <- intraspecies_synteny(blast_list, pdata$annotation,anchors=anchors, max_gaps=max_gaps, e_value=value, intra_dir = "data_nco/nco_25_5_1e-5/" )

nco_syn <- parse_collinearity(intrasyn, as = "all")
write.table(nco_syn, file='data_nco/nco_25_5_1e-5/intra_25_5_1e-5.tsv', quote=FALSE, sep='\t', row.names = FALSE)




value = 1e-5
anchors = 5
max_gaps = 35

intrasyn <- intraspecies_synteny(blast_list, pdata$annotation,anchors=anchors, max_gaps=max_gaps, e_value=value, intra_dir = "data_nco/nco_35_5_1e-5/" )

nco_syn <- parse_collinearity(intrasyn, as = "all")
write.table(nco_syn, file='data_nco/nco_35_5_1e-5/intra_35_5_1e-5.tsv', quote=FALSE, sep='\t', row.names = FALSE)





value = 1e-5
anchors = 3
max_gaps = 35

intrasyn <- intraspecies_synteny(blast_list, pdata$annotation,anchors=anchors, max_gaps=max_gaps, e_value=value, intra_dir = "data_nco/nco_35_3_1e-5/" )

nco_syn <- parse_collinearity(intrasyn, as = "all")
write.table(nco_syn, file='data_nco/nco_35_3_1e-5/intra_35_3_1e-5.tsv', quote=FALSE, sep='\t', row.names = FALSE)


value = 1e-5
anchors = 3
max_gaps = 50

intrasyn <- intraspecies_synteny(blast_list, pdata$annotation,anchors=anchors, max_gaps=max_gaps, e_value=value, intra_dir = "data_nco/nco_50_3_1e-5/" )

nco_syn <- parse_collinearity(intrasyn, as = "all")
write.table(nco_syn, file='data_nco/nco_50_3_1e-5/intra_50_3_1e-5.tsv', quote=FALSE, sep='\t', row.names = FALSE)


value = 1e-3
anchors = 3
max_gaps = 35

intrasyn <- intraspecies_synteny(blast_list, pdata$annotation,anchors=anchors, max_gaps=max_gaps, e_value=value, intra_dir = "data_nco/nco_35_3_1e-3/" )

nco_syn <- parse_collinearity(intrasyn, as = "all")
write.table(nco_syn, file='data_nco/nco_35_3_1e-3/intra_35_3_1e-3.tsv', quote=FALSE, sep='\t', row.names = FALSE)


value = 1e-3
anchors = 3
max_gaps = 55

intrasyn <- intraspecies_synteny(blast_list, pdata$annotation,anchors=anchors, max_gaps=max_gaps, e_value=value, intra_dir = "data_nco/nco_55_3_1e-3/" )

nco_syn <- parse_collinearity(intrasyn, as = "all")
write.table(nco_syn, file='data_nco/nco_55_3_1e-3/intra_55_3_1e-3.tsv', quote=FALSE, sep='\t', row.names = FALSE)




####Inter synteny

##ambo ive

library(syntenet)
library(GenomicRanges)

Sys.setenv(
    PATH = paste(
        Sys.getenv("PATH"), "/beegfs/data/mmorel/miniconda3/bin", sep = ":"
    )
)
setwd("/beegfs/data/mmorel/Infer_Anth/data/")

aastringsetlist <- fasta2AAStringSetlist("inter_species/")

annotation_illicium = read.csv("inter_species/illicium_verum_no_splicing.tsv", sep="\t")
annotation_amborella = read.csv("inter_species/amborella_trichopoda_nosplicing_63666.tsv", sep="\t")

annot_ive = makeGRangesFromDataFrame(annotation_illicium, seqnames.field="seqnames", start.field = "start", 
                                 end.field="end", strand.field="strand", keep.extra.columns=TRUE)

annot_ambo = makeGRangesFromDataFrame(annotation_amborella, seqnames.field="seqnames", start.field = "start", 
                                 end.field="end", strand.field="strand", keep.extra.columns=TRUE)


annot_lst= list(annot_ambo, annot_ive) #ordre alphabétique des diminutifs
names(annot_lst) <- c("atr", "ive") #ordre alphabétique
names(aastringsetlist) <- c("atr", "ive")
check_input(aastringsetlist, annot_lst)


#########

pdata <- process_input(aastringsetlist, annot_lst)

blast_list <- run_diamond(seq = pdata$seq, top_hits=5)
write.table(do.call(rbind, blast_list), file='inter_species/ambo_illicium/blast_table_atr_ive.tsv', quote=FALSE, sep='\t', row.names = FALSE)

######


diamond_inter <- blast_list[c(2,3,4,7,8,12)]

value = 1e-3
anchors = 3
max_gaps = 55


intersyn <- interspecies_synteny(diamond_inter, pdata$annotation, inter_dir = "inter_species/ambo_illicium/55_3_1e-3/" ,anchors=anchors, max_gaps=max_gaps, e_value=value)
#.collinearity

all_syn <- parse_collinearity(intersyn, as = "all")
write.table(all_syn , file='inter_species/ambo_illicium/55_3_1e-3/inter_55_3_1e-3.tsv', quote=FALSE, sep='\t', row.names = FALSE)







####Inter synteny


##aristo vs ambo vs vitis vs nymphaea vs aco vs ive
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
#annotation_cycas = read.csv("inter_species/Cycas_panzhihuaensis.tsv", sep="\t")
#annotation_annona = read.csv("inter_species/annona_cherimola_no_splicing.tsv", sep="\t")
annotation_contorta = read.csv("inter_species/aristo_contorta_no_splicing.tsv", sep="\t")
annotation_nymphaea = read.csv("inter_species/nymphaea_colorata.tsv", sep="\t")
annotation_illicium = read.csv("inter_species/illicium_verum_no_splicing.tsv", sep="\t")

annot_ambo = makeGRangesFromDataFrame(annotation_ambo, seqnames.field="seqnames", start.field = "start", 
                         end.field="end", strand.field="strand", keep.extra.columns=TRUE)
annot_aristo = makeGRangesFromDataFrame(annotation_aristo, seqnames.field="seqnames", start.field = "start", 
                                 end.field="end", strand.field="strand", keep.extra.columns=TRUE)
annot_vitis = makeGRangesFromDataFrame(annotation_vitis, seqnames.field="seqnames", start.field = "start", 
                                 end.field="end", strand.field="strand", keep.extra.columns=TRUE)
#annot_cycas = makeGRangesFromDataFrame(annotation_cycas, seqnames.field="seqnames", start.field = "start", 
#                                 end.field="end", strand.field="strand", keep.extra.columns=TRUE)
#annot_annona = makeGRangesFromDataFrame(annotation_annona, seqnames.field="seqnames", start.field = "start", 
#                                 end.field="end", strand.field="strand", keep.extra.columns=TRUE)
annot_contorta = makeGRangesFromDataFrame(annotation_contorta, seqnames.field="seqnames", start.field = "start", 
                                 end.field="end", strand.field="strand", keep.extra.columns=TRUE)

annot_nymphaea = makeGRangesFromDataFrame(annotation_nymphaea, seqnames.field="seqnames", start.field = "start", 
                                 end.field="end", strand.field="strand", keep.extra.columns=TRUE)
annot_illicium = makeGRangesFromDataFrame(annotation_illicium, seqnames.field="seqnames", start.field = "start", 
                                 end.field="end", strand.field="strand", keep.extra.columns=TRUE)

#annot_lst= list(annot_annona, annot_contorta, annot_aristo, annot_ambo, annot_cycas, annot_vitis) #ordre alphabétique des diminutifs
annot_lst= list(annot_contorta, annot_aristo, annot_ambo, annot_illicium, annot_nymphaea, annot_vitis) #ordre alphabétique des diminutifs
#names(annot_lst) <- c("ach", "aco", "afi", "atr", "cyp", "viv") #ordre alphabétique
names(annot_lst) <- c("aco","afi", "atr","ive", "nco", "viv") #ordre alphabétique
#names(aastringsetlist) <- c("ach", "aco", "afi", "atr", "cyp", "viv")
names(aastringsetlist) <- c("aco", "afi", "atr","ive", "nco", "viv")
check_input(aastringsetlist, annot_lst)


#########

pdata <- process_input(aastringsetlist, annot_lst)

#blast_list <- run_diamond(seq = pdata$seq, top_hits=5) #, evalue=0.001
#write.table(do.call(rbind, blast_list), file='inter_species/aco_afi_atr_ive_nco_viv/blast_table.tsv', quote=FALSE, sep='\t', row.names = FALSE)

######
#to not have to run it multiple times
blast_list = read.table('inter_species/aco_afi_atr_ive_nco_viv/blast_table.tsv', sep = '\t',quote=FALSE)



'''

> names(blast_list)
 [1] "aco_aco" "aco_afi" "aco_atr" "aco_ive" "aco_nco" "aco_viv" "afi_aco"
 [8] "afi_afi" "afi_atr" "afi_ive" "afi_nco" "afi_viv" "atr_aco" "atr_afi"
[15] "atr_atr" "atr_ive" "atr_nco" "atr_viv" "ive_aco" "ive_afi" "ive_atr"
[22] "ive_ive" "ive_nco" "ive_viv" "nco_aco" "nco_afi" "nco_atr" "nco_ive"
[29] "nco_nco" "nco_viv" "viv_aco" "viv_afi" "viv_atr" "viv_ive" "viv_nco"
[36] "viv_viv"

'''



'''
> names(blast_list)
 [1] "aco_aco" "aco_afi" "aco_atr" "aco_nco" "aco_viv" "afi_aco" "afi_afi"
 [8] "afi_atr" "afi_nco" "afi_viv" "atr_aco" "atr_afi" "atr_atr" "atr_nco"
[15] "atr_viv" "nco_aco" "nco_afi" "nco_atr" "nco_nco" "nco_viv" "viv_aco"
[22] "viv_afi" "viv_atr" "viv_nco" "viv_viv"
'''

diamond_inter <- blast_list[c(2,3,4,5,8,9,10,14,15,20)]


"aco_afi" "aco_atr" "aco_nco" "aco_viv"
"afi_atr" "afi_nco" "afi_viv"
"atr_nco" "atr_viv" 
"nco_viv"


#diamond_inter <- blast_list[c(2,3,4,5,6,9,10,11,12,16,17,18,23,24,30)]

diamond_inter <- blast_list[c(2,3,4,7,8,12)]

value = 1e-5
anchors = 5
max_gaps = 25

intersyn <- interspecies_synteny(diamond_inter, pdata$annotation, inter_dir = "inter_species/aco_afi_atr_ive_nco_viv/25_5_1e-5/" ,anchors=anchors, max_gaps=max_gaps, e_value=value)
#.collinearity

all_syn <- parse_collinearity(intersyn, as = "all")
write.table(all_syn , file='inter_species/aco_afi_atr_ive_nco_viv/25_5_1e-5/inter_25_5_1e-5.tsv', quote=FALSE, sep='\t', row.names = FALSE)



# value = 1e-5
# anchors = 5
# max_gaps = 25

# intersyn <- interspecies_synteny(diamond_inter, pdata$annotation, inter_dir = "inter_species/ambo_aristo_nymphaea_vitis/25_5_1e-5/" ,anchors=anchors, max_gaps=max_gaps, e_value=value)
# #.collinearity

# all_syn <- parse_collinearity(intersyn, as = "all")
# write.table(all_syn , file='inter_species/ambo_aristo_nymphaea_vitis/25_5_1e-5/inter_25_5_1e-5.tsv', quote=FALSE, sep='\t', row.names = FALSE)



# value = 1e-5
# anchors = 5
# max_gaps = 25

# intersyn <- interspecies_synteny(diamond_inter, pdata$annotation, inter_dir = "inter_species/ambo_aristo_nymphaea_vitis/25_5_1e-5/" ,anchors=anchors, max_gaps=max_gaps, e_value=value)
# #.collinearity

# all_syn <- parse_collinearity(intersyn, as = "all")
# write.table(all_syn , file='inter_species/ambo_aristo_nymphaea_vitis/25_5_1e-5/inter_25_5_1e-5.tsv', quote=FALSE, sep='\t', row.names = FALSE)




value = 1e-5
anchors = 5
max_gaps = 35

intersyn <- interspecies_synteny(diamond_inter, pdata$annotation, inter_dir = "inter_species/ambo_aristo_nymphaea_vitis/35_5_1e-5/" ,anchors=anchors, max_gaps=max_gaps, e_value=value)
#.collinearity

all_syn <- parse_collinearity(intersyn, as = "all")
write.table(all_syn , file='inter_species/ambo_aristo_nymphaea_vitis/35_5_1e-5/inter_35_5_1e-5.tsv', quote=FALSE, sep='\t', row.names = FALSE)



value = 1e-5
anchors = 3
max_gaps = 35

intersyn <- interspecies_synteny(diamond_inter, pdata$annotation, inter_dir = "inter_species/ambo_aristo_nymphaea_vitis/35_3_1e-5/" ,anchors=anchors, max_gaps=max_gaps, e_value=value)
#.collinearity

all_syn <- parse_collinearity(intersyn, as = "all")
write.table(all_syn , file='inter_species/ambo_aristo_nymphaea_vitis/35_3_1e-5/inter_35_3_1e-5.tsv', quote=FALSE, sep='\t', row.names = FALSE)


value = 1e-5
anchors = 3
max_gaps = 50


intersyn <- interspecies_synteny(diamond_inter, pdata$annotation, inter_dir = "inter_species/ambo_aristo_nymphaea_vitis/50_3_1e-5/" ,anchors=anchors, max_gaps=max_gaps, e_value=value)
#.collinearity

all_syn <- parse_collinearity(intersyn, as = "all")
write.table(all_syn , file='inter_species/ambo_aristo_nymphaea_vitis/50_3_1e-5/inter_50_3_1e-5.tsv', quote=FALSE, sep='\t', row.names = FALSE)


value = 1e-3
anchors = 3
max_gaps = 35


intersyn <- interspecies_synteny(diamond_inter, pdata$annotation, inter_dir = "inter_species/ambo_aristo_nymphaea_vitis/35_3_1e-3/" ,anchors=anchors, max_gaps=max_gaps, e_value=value)
#.collinearity

all_syn <- parse_collinearity(intersyn, as = "all")
write.table(all_syn , file='inter_species/ambo_aristo_nymphaea_vitis/35_3_1e-3/inter_35_3_1e-3.tsv', quote=FALSE, sep='\t', row.names = FALSE)


value = 1e-3
anchors = 3
max_gaps = 55


intersyn <- interspecies_synteny(diamond_inter, pdata$annotation, inter_dir = "inter_species/aco_afi_atr_ive_nco_viv/55_3_1e-3/" ,anchors=anchors, max_gaps=max_gaps, e_value=value)
#.collinearity

all_syn <- parse_collinearity(intersyn, as = "all")
write.table(all_syn , file='inter_species/aco_afi_atr_ive_nco_viv/55_3_1e-3/inter_55_3_1e-3.tsv', quote=FALSE, sep='\t', row.names = FALSE)



yabby_afi = list("afi_afi_Af03G134500", "afi_afi_Af01G026600", "afi_afi_Af02G072800", "afi_afi_Af05G185300","afi_afi_Af07G024800")
#INO, YAB5, FIL, YAB2, CRC
yabby_atr = list("atr_atr_94506842", "atr_atr_94491478", "atr_atr_94503013", "atr_atr_94473924", "atr_atr_94507018")
#YAB2, CRC, INO, YAB5, FIL
yabby_viv = list("viv_Vitvi06g00972", "viv_Vitvi02g00510", "viv_Vitvi15g00708", "viv_Vitvi11g00492", "viv_Vitvi01g00013",  "viv_Vitvi08g00274", "viv_viv_Vitvi01g00703")
# YAB2 FIL FIL YAB5 CRC YAB2 INO
yabby_cycas = list("cyp_CYCAS_007102", "cyp_CYCAS_008411", "cyp_CYCAS_022386")


afi_afi_Af03G134500 afi_afi_Af01G026600 afi_afi_Af02G072800 afi_afi_Af05G185300 afi_afi_Af07G024800

nco_Nycol.B00216.v1 nco_Nycol.N00117.v1 nco_Nycol.L00815.v1 nco_Nycol.D01480.v1 nco_Nycol.A00487.v1 nco_Nycol.G00354.v1
nco_GWHGAAYW018449 nco_GWHGAAYW010642 nco_GWHGAAYW008701 nco_GWHGAAYW012319 nco_GWHGAAYW011460 nco_GWHGAAYW000547 nco_GWHGAAYW023292




####Inter synteny


##aristo contorta vs nymphaea
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

annotation_contorta = read.csv("inter_species/aristo_contorta_no_splicing.tsv", sep="\t")
annotation_nymphaea = read.csv("inter_species/nymphaea_colorata.tsv", sep="\t")

annot_contorta = makeGRangesFromDataFrame(annotation_contorta, seqnames.field="seqnames", start.field = "start", 
                                 end.field="end", strand.field="strand", keep.extra.columns=TRUE)

annot_nymphaea = makeGRangesFromDataFrame(annotation_nymphaea, seqnames.field="seqnames", start.field = "start", 
                                 end.field="end", strand.field="strand", keep.extra.columns=TRUE)


annot_lst= list(annot_contorta, annot_nymphaea) #ordre alphabétique des diminutifs
names(annot_lst) <- c("aco", "nco") #ordre alphabétique
names(aastringsetlist) <- c("aco", "nco")
check_input(aastringsetlist, annot_lst)


#########

pdata <- process_input(aastringsetlist, annot_lst)

blast_list <- run_diamond(seq = pdata$seq, top_hits=5)
write.table(do.call(rbind, blast_list), file='inter_species/contorta_nymphaea/blast_table_aco_nco.tsv', quote=FALSE, sep='\t', row.names = FALSE)

######


diamond_inter <- blast_list[c(2,3,4,7,8,12)]

value = 1e-3
anchors = 3
max_gaps = 55


intersyn <- interspecies_synteny(diamond_inter, pdata$annotation, inter_dir = "inter_species/contorta_nymphaea/55_3_1e-3/" ,anchors=anchors, max_gaps=max_gaps, e_value=value)
#.collinearity

all_syn <- parse_collinearity(intersyn, as = "all")
write.table(all_syn , file='inter_species/contorta_nymphaea/55_3_1e-3/inter_55_3_1e-3.tsv', quote=FALSE, sep='\t', row.names = FALSE)




# value = 1e-3
# anchors = 3
# max_gaps = 55
# 
# intersyn <- interspecies_synteny(diamond_inter, pdata$annotation, inter_dir = "inter_species/ambo_aristo_nymphaea_vitis/25_5_1e-5/" ,anchors=anchors, max_gaps=max_gaps, e_value=value)
#.collinearity
# 
# all_syn <- parse_collinearity(intersyn, as = "all")
# write.table(all_syn , file='inter_species/ambo_aristo_nymphaea_vitis/25_5_1e-5/inter_25_5_1e-5.tsv', quote=FALSE, sep='\t', row.names = FALSE)

####Inter synteny


##ambo vs vitis

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

annotation_contorta = read.csv("inter_species/aristo_contorta_no_splicing.tsv", sep="\t")
annotation_nymphaea = read.csv("inter_species/nymphaea_colorata.tsv", sep="\t")

annot_contorta = makeGRangesFromDataFrame(annotation_contorta, seqnames.field="seqnames", start.field = "start", 
                                 end.field="end", strand.field="strand", keep.extra.columns=TRUE)

annot_nymphaea = makeGRangesFromDataFrame(annotation_nymphaea, seqnames.field="seqnames", start.field = "start", 
                                 end.field="end", strand.field="strand", keep.extra.columns=TRUE)


annot_lst= list(annot_contorta, annot_nymphaea) #ordre alphabétique des diminutifs
names(annot_lst) <- c("aco", "nco") #ordre alphabétique
names(aastringsetlist) <- c("aco", "nco")
check_input(aastringsetlist, annot_lst)


#########

pdata <- process_input(aastringsetlist, annot_lst)

blast_list <- run_diamond(seq = pdata$seq, top_hits=5)
write.table(do.call(rbind, blast_list), file='inter_species/contorta_nymphaea/blast_table_aco_nco.tsv', quote=FALSE, sep='\t', row.names = FALSE)

######


diamond_inter <- blast_list[c(2,3,4,7,8,12)]

value = 1e-3
anchors = 3
max_gaps = 55


intersyn <- interspecies_synteny(diamond_inter, pdata$annotation, inter_dir = "inter_species/contorta_nymphaea/55_3_1e-3/" ,anchors=anchors, max_gaps=max_gaps, e_value=value)
#.collinearity

all_syn <- parse_collinearity(intersyn, as = "all")
write.table(all_syn , file='inter_species/contorta_nymphaea/55_3_1e-3/inter_55_3_1e-3.tsv', quote=FALSE, sep='\t', row.names = FALSE)