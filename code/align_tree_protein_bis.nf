//path_file = "/beegfs/data/mmorel/Infer_Anth/datation/"
//path_file = "/beegfs/data/mmorel/catarina_gls/update_dataset/"
path_file = "/beegfs/data/mmorel/Infer_Anth/align/"

params.sequences = path_file + "gymno_nooutgroup.nt.fa"
//"GLS_cds.rename.fa"
//seed_plants
//"cds_only_catalytic_domain.fst"
//params.topology = path_file + "test_constrained_salse_nopara.nwk"
params.resdir= path_file +  "gymno/no_outgroup/results_gymno_aa/" 
//"results_outgroup/" //do not forget to change the resdir
//params.clean = "true"// trimal, keep_all, true. 
//params.clean_aa = "hmm" //"hmm", "goalign"
params.clean_trimal = "gappy_02" //gappy_09 (=sites with 10% gaps), auto, gappy_02 (=sites with 80% gaps) 
params.tree = "true" //to construct the tree at the end
params.hom_macse = "false" //if we trim non homologuous fragment at the beginning
params.guide_tree = "false" //if we want to constrain the topology true|false
params.complextree = "false"

sequences= file(params.sequences)
//topology = file(params.topology)
//clean = params.clean
//clean_aa = params.clean_aa
clean_trimal = params.clean_trimal
tree = params.tree
guide_tree = params.guide_tree
hom_macse = params.hom_macse
complextree = params.complextree
resdir=file(params.resdir)
resdir.with {mkdirs()}
// récupérer toutes les séquences du fichier donné par Charlie

//trim non homologuous sequences
process trim_nonhomologous{
    publishDir "${resdir}", mode: 'copy'
    label 'tool'
    input: 
    file sequences
    output:
    file "*trim*.nt.fasta" into TrimCoding
    file "*trim*.aa.fasta" into TrimTranslated

    shell:
    if ( hom_macse == "true")
    '''
    java -jar ~/bin/macse_v2.07.jar -prog trimNonHomologousFragments -seq !{sequences} -out_trim_info output_stats.csv \
    -out_NT trim_yabby.nt.fasta -out_AA trim_yabby.aa.fasta
    '''  
    else if (hom_macse == "false")
    '''
    goalign translate --unaligned -i !{sequences} -o notrim_yabby.aa.fasta
    cp  !{sequences} notrim_yabby.nt.fasta
    '''
}

//align homologuous sequences with various aligners
process align_macse{
    label 'tool'
    //publishDir "${resdir}", mode: 'copy'
    input:
    file trim_sequences from TrimCoding
    output:
    file 'macse.yabby_sequences.aa.fasta' into MACSE_Align
    shell:
    //aligns nucleotide (NT) coding sequences using their amino acid (AA) translations.
    '''
    java -jar ~/bin/macse_v2.07.jar -prog alignSequences -seq !{trim_sequences} -out_NT macse.yabby_sequences.fasta -out_AA macse.yabby_sequences.aa.fasta
    '''
}


//protein alignment using mafft
process align_mafft{
    label 'tool'

    //publishDir "${resdir}", mode: 'copy'
    input:
    file seq from TrimTranslated
    output:
    file 'mafft.*' into MafftAlign
    shell:
    '''
    mafft --maxiterate 1000 --localpair !{seq} > mafft.!{seq}
    '''
}

//protein alignment using muscle
process align_muscle{
    label 'tool'
    //publishDir "${resdir}", mode: 'copy'
    input:
    file seq from TrimTranslated
    output:
    file('muscle.*') into MuscleAlign
    shell:
    '''
    muscle -align !{seq} -output muscle.!{seq}
    '''
}

//protein alignment using clustalo
process align_clustal{
    label 'tool'
    //publishDir "${resdir}", mode: 'copy'
    input:
    file seq from TrimTranslated
    output:
    file 'clustal*.fasta' into ClustalAlign
    shell:
    '''
    clustalo -i !{seq} -o clustal.!{seq} --iter=2
    '''
}

// //protein alignment using prank
// process align_prank{
//     //publishDir "${resdir}", mode: 'copy'
//     input:
//     file seq from TrimTranslated
//     output:
//     file 'prank*.fasta' into PrankAlign
//     shell:
//     '''
//     prank -d=!{seq} -o=prank.!{seq} -F
//     mv prank.!{seq}.best.fas prank.!{seq}
//     '''
// }
//PrankAlign

AllAlignAA = ClustalAlign.mix(MACSE_Align,MuscleAlign,MafftAlign)

AllAlignAA.into {
AllAlign1
AllAlign2
AllAlign3
}


process clean_hmm{
    label 'tool'
    publishDir "${resdir}", mode: 'copy'
    input: 
    file(align_aa) from AllAlign1
    output:
    file("*hmm.fasta") into HmmClean_AA
    //when : clean == "true"
    shell:
    '''
    apptainer exec --bind /beegfs --home $HOME:/home/$USER /beegfs/home/mmorel/.Singularity_NF/hmmcleaner.sif HmmCleaner.pl !{align_aa} --large
    '''
}

process clean_trimal{

    //publishDir "${resdir}", mode: 'copy'
    input: 
    file (align_aa) from AllAlign2
    output:
    file("trimal_clean.*") into Trimalclean_AA //trim align

    //when : clean == "trimal" 

    shell:
    if( clean_trimal == 'auto' )
    '''
    trimal -in !{align_aa} -out trimal_clean.!{align_aa} -automated1
    '''
    else if( clean_trimal == 'gappy_09' ) //will remove all columns with gaps in more than 10% of the sequences 
    '''
    trimal -in !{align_aa} -out trimal_clean.!{align_aa} -gt 0.9  
    '''
    else if( clean_trimal == 'gappy_02' ) //will remove all columns with gaps in more than 80% of the sequences 
    '''
    trimal -in !{align_aa} -out trimal_clean.!{align_aa} -gt 0.2  
    '''
    else 
    error "Invalid alignment mode: ${clean_trimal}"
}

//trimal -in !{align} -out trimal.!{align} -automated1
//trimal -in !{align} -backtrans !{dna_seq} -out clean_codon.!{align} -automated1
//trimal -in !{align} -backtrans !{dna_seq} -out clean_codon.!{align} -gt 0.9
//trimal -in !{align} -backtrans !{dna_seq} -out clean_codon.!{align} -selectcols {0-88, 155-191, 255-323}



Clean_or_not = AllAlign3.mix(Trimalclean_AA, HmmClean_AA)


process reconstruct_tree{
    label 'iqtree2'
    publishDir "${resdir}", mode: 'copy'
    input:
    file(align_aa) from Clean_or_not
    //file topology
    output:
    tuple file(align_aa), file('*.treefile') into TreeFile

    when : tree == "true" 
    shell:
    if (guide_tree == "false")
    '''
    iqtree -s !{align_aa} -m LG+R7+F -nt !{task.cpus} -B 1000 -safe -bnni
    '''
    else if (guide_tree == "true")
    '''
    iqtree -s !{align_aa} -m LG+R7 -g !{topology} -nt !{task.cpus} -B 1000 -safe -bnni
    '''
}

process reconstruct_tree_guided{
    label 'iqtree2'
    publishDir "${resdir}", mode: 'copy'
    input:
    tuple file(align_aa), file(guidetree) from TreeFile
    //file topology
    output:
    tuple file(align_aa), file('complex*.treefile'), file('complex*.iqtree')

    when : complextree == "true" 
    shell:
    '''
    s=!{align_aa}
    name="${s%.*}"
    iqtree -s !{align_aa} -m LG+C60+R4 -ft !{guidetree} -nt !{task.cpus} --prefix complex_${name}
    '''
    
}
