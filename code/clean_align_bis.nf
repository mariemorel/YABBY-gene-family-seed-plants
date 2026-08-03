
//path_file = "/beegfs/data/mmorel/Infer_Anth/datation/"
path_file = "/beegfs/data/mmorel/Infer_Anth/align/"

params.sequences = path_file + "seed_plants_cycas_domains_nt.fa"
//params.sequences = path_file + "samples1_seqs_dna_gymno_noB.fasta"
//params.topology = path_file + "test_constrained_salse_nopara.nwk"
params.resdir= path_file + "seedplants/cycas/results_seedplants_cycas_nt_domain/" //do not forget to change the resdir

//params.clean = "true"// trimal, keep_all, true. 
//params.clean_aa = "hmm" //"hmm", "goalign"
params.clean_trimal = "gappy_02" //gappy_09 (=sites with 10% gaps), auto, gappy_01 (=sites with 90% gaps) 
params.tree = "true" //to construct the tree at the end
params.hom_macse = "false" //if we trim non homologuous fragment at the beginning
params.guide_tree = "false" //if we want to constrain the topology true|false
params.outgroup = path_file + "outgroup.txt"

sequences= file(params.sequences)
outgroup = file(params.outgroup)
//topology = file(params.topology)
//clean = params.clean
//clean_aa = params.clean_aa
clean_trimal = params.clean_trimal
tree = params.tree
guide_tree = params.guide_tree
hom_macse = params.hom_macse
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
    file "*trim*.nt.fasta" into TrimCoding, TrimCoding1, TrimCoding2
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
    tuple file('codon.macse.yabby_sequences.aa.fasta'), file('codon.macse.yabby_sequences.fasta') into MACSE_Align
    shell:
    '''
    java -jar ~/bin/macse_v2.07.jar -prog alignSequences -seq !{trim_sequences} -out_NT codon.macse.yabby_sequences.fasta -out_AA codon.macse.yabby_sequences.aa.fasta
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

AllAlignAA = ClustalAlign.mix(MuscleAlign,MafftAlign)

process codon_transfo{
    label 'tool'

    //publishDir "${resdir}", mode: 'copy'
    input: 
    file(align_nt) from TrimCoding1
    file(align_aa) from AllAlignAA
    
    output:
    tuple file(align_aa), file("codon*") into TransfoCodon //align in codon after aa clean

    shell:
    '''
    goalign codonalign -f !{align_nt} -i !{align_aa} -o codon.!{align_aa}
    '''
}


AllAlign = TransfoCodon.mix(MACSE_Align)
AllAlign.into {
AllAlign1
AllAlign2
AllAlign3
}



process clean_hmm{
    label 'tool'
    //publishDir "${resdir}", mode: 'copy'
    input: 
    tuple file(align_aa), file(align_codon) from AllAlign1
    output:
    tuple file("*hmm.fasta"), file(align_codon) into Clean_AA
    //when : clean == "true"
    shell:
    '''
    apptainer exec --bind /beegfs --home $HOME:/home/$USER /beegfs/home/mmorel/.Singularity_NF/hmmcleaner.sif HmmCleaner.pl !{align_aa}
    '''
}

process codon_mask{
    //publishDir "${resdir}", mode: 'copy'
    label 'tool'
    input: 
    tuple file(align_aa), file(align_codon) from Clean_AA
    output:
    tuple file(align_aa), file("transfo*.fasta") into CodonAlign
  
    shell:
    '''
    java -jar ~/bin/macse_v2.07.jar -prog reportMaskAA2NT -align !{align_codon} -align_AA !{align_aa} -mask_AA - -out_NT transfo_!{align_codon} -out_mask_detail mask_detail_!{align_codon}
    '''  
}



//also possible to trim in codon directly with trimal

process clean_trimal{

    //publishDir "${resdir}", mode: 'copy'
    input: 
    file dna_seq from TrimCoding2
    tuple file(align_aa), file(align_codon) from AllAlign2
    output:
    tuple file(align_aa), file("trimal_clean.*") into TrimalCodon //align in codon + trim

    //when : clean == "trimal" 

    shell:
    if( clean_trimal == 'auto' )
    '''
    trimal -in !{align_aa} -backtrans !{dna_seq} -out trimal_clean.!{align_codon} -automated1
    '''
    else if( clean_trimal == 'gappy_09' ) //will remove all columns with gaps in more than 10% of the sequences 
    '''
    trimal -in !{align_aa} -backtrans !{dna_seq} -out trimal_clean.!{align_codon} -gt 0.9  
    '''
    else if( clean_trimal == 'gappy_02' ) //will remove all columns with gaps in more than 90% of the sequences 
    '''
    trimal -in !{align_aa} -backtrans !{dna_seq} -out trimal_clean.!{align_codon} -gt 0.2  
    '''
    else 
    error "Invalid alignment mode: ${clean_trimal}"
}

//trimal -in !{align} -out trimal.!{align} -automated1
//trimal -in !{align} -backtrans !{dna_seq} -out clean_codon.!{align} -automated1
//trimal -in !{align} -backtrans !{dna_seq} -out clean_codon.!{align} -gt 0.9
//trimal -in !{align} -backtrans !{dna_seq} -out clean_codon.!{align} -selectcols {0-88, 155-191, 255-323}



//retrieve the align in codon either after transfo or not

Clean_all = AllAlign3.mix(TrimalCodon, CodonAlign)

process reconstruct_tree{
    label 'iqtree2'
    publishDir "${resdir}", mode: 'copy'
    input:
    tuple file(align_aa), file(align_codon) from Clean_all
    //file topology
    output:
    tuple file(align_codon), file('*.treefile'), file('*.iqtree') into Consensus

    when : tree == "true" 
    shell:
    if (guide_tree == "false")
    '''
    length=$(awk '/^>/ {if(found) exit; found=1; next} found {len += length} END {print len}' !{align_codon})

    echo "#nexus  
    BEGIN SETS;
    charset !{align_codon} = 1-${length}\\3, 2-${length}\\3; 
    END;" > partitions.nx
    
    iqtree -s !{align_codon} -st DNA -p partitions.nx -m GTR+R4+F+I -nt !{task.cpus} -B 1000 -safe -bnni --prefix !{align_codon}
    '''
    else if (guide_tree == "true")
    '''
    length=$(awk '/^>/ {if(found) exit; found=1; next} found {len += length} END {print len}' !{align_codon})
    
    echo "#nexus  
    BEGIN SETS;
    charset !{align_codon} = 1-${length}\\3, 2-${length}\\3; 
    END;" > partitions.nx

    iqtree -s !{align_codon} -st DNA -g !{topology} -p partitions.nx -m GTR+R4+F+I -nt !{task.cpus} -B 1000 -safe -bnni --prefix !{align_codon}
    '''
}


AllTrees = Consensus.flatten().collect()

process consensus_tree{
    label 'tool'
    publishDir "${resdir}", mode: 'copy'
    input:
    file("*") from AllTrees
    file outgroup
    output:
    file("*.consensus")

    //for i in `ls transfo_codon*.treefile.reroot` ; do cat $i >> reroot.hmm_tree.tree ; done
    '''
    for i in `ls *.treefile` ; do gotree reroot outgroup -l !{outgroup} -i $i > reroot.$i ; done
    for i in `ls codon*.treefile.reroot` ; do cat $i >> reroot.notrim_tree.tree ; done
    for i in `ls trimal_clean*.treefile.reroot` ; do cat $i >> reroot.trimal_tree.tree ; done
    for i in reroot.notrim_tree.tree reroot.trimal_tree.tree ; do gotree compute consensus -f 0.5 -i $i > $i.consensus; done
    '''

}


/*  

//baliphy
//faire en sorte de le lancer plusieurs fois each x from

process baliphy{
    conda '/beegfs/data/mmorel/miniconda3/envs/baliphy'
    publishDir "${resdir}", mode: 'copy'
    input:
    each x from {1..4}
    file sequences
    output:
    
    shell:
    '''
    bali-phy !{sequences} -S gy94[pi=f1x4] -A Codons 
    '''
}

*/


/*
process clean_hmmerclean{
    
    publishDir "${resdir}", mode: 'copy'
    input: 
    file clean_align from AllAlign1
    output:
    file "*hmm.fasta" into CleanAlignAA0 //align in aa

    when : clean == "hmm"
    shell:
    '''
    apptainer exec --bind /beegfs --home $HOME:/home/$USER /beegfs/home/mmorel/.Singularity_NF/hmmcleaner.sif HmmCleaner.pl !{clean_align}
    '''
}

process clean_goalign{
    label 'goalign'

    publishDir "${resdir}", mode: 'copy'
    input: 
    file align from AllAlign2
    output:
    file "clean.*" into CleanAlignAA //align in aa
  
    when : clean == "goalign" 
   
    shell:
    '''
    goalign clean sites -c 0.1 -i !{align} -o clean.!{align}
    '''
}

//Clean_AA = CleanAlignAA0.concat(CleanAlignAA) //retrieve either goalign clean or hmm cleaner. Need to translate it back into codon.
*/
