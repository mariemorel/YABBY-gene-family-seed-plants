path_file = "/beegfs/data/mmorel/Infer_Anth/align/seedplants/cycas/results_seedplants_cycas_codon_domain/" 

params.resdir = path_file + "annotated_trees/"
params.outgroup = "/beegfs/data/mmorel/Infer_Anth/align/outgroup.txt" 

outgroup = file(params.outgroup)

resdir=file(params.resdir)
resdir.with {mkdirs()}

trees = Channel.fromPath(path_file + '*.treefile').set{Trees_ch}


process collapse_tree{
    label 'tool'
    input: 
    file outgroup
    file tree from Trees_ch
    output:
    set file(tree), file ("collapse*.treefile") into AllTrees
    shell:
    '''
    gotree reroot outgroup -i !{tree} -l !{outgroup} > root.!{tree}
    gotree collapse support -s 50 -i root.!{tree} > collapse.!{tree}
    '''
}



process comment_bootstrap{
    label 'tool'
    input: 
    file tree from AllTrees.flatten()
    output:
    file "rename.*.tree" into Prepared_Tree
    shell:
    '''
    NAME=!{tree}
    core_name="${NAME%.*}"
    if grep -E ')[0-9]{1,3}[:]' !{tree} > /dev/null; 
    then
        cat !{tree} | sed -E 's/\\)([0-9]+):/)"\\1":/g' > quoted_bootstrap.tree ;
        gotree comment transfer -i quoted_bootstrap.tree > commented_bootstrap.tree ;
        gotree rename --internal --auto --tips=false -i commented_bootstrap.tree > rename.${core_name}.tree; 
    else gotree rename --internal --auto --tips=false -i !{tree} > rename.${core_name}.tree; 
    fi   
    '''
}


process annotation_file{
    label 'python'
    input:
    file tree from Prepared_Tree
    output:
    set file(tree), file("internal_nodes*") into Color_file
    shell:
    '''
    annotate_yabby.py !{tree}
    '''
}


process color_tree{
    label 'tool'
    publishDir "${resdir}", mode: 'copy'

    input:
    set file(tree), file(annot) from Color_file
    output:
    file '*.nx'

    shell:
    '''
    gotree comment clear -i !{tree} | gotree rename -m !{annot} > labeled.tree
    #remove internal node names
    sed -e 's/N00000[0-9]*//g' labeled.tree |  gotree reformat nexus > !{tree}.nx
    '''
}