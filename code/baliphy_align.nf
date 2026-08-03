path_file = "/beegfs/data/mmorel/Infer_Anth/align/"

params.sequences = path_file + "gymno_dataset_nt.cycas.fa"
//"subset_seed_plants.nt.fa"
// "seqs_dna_update_nyth_less_seq.fa"
//params.resdir= path_file + "results_baliphy_codon_outgroup/" //do not forget to change the resdir

sequences= file(params.sequences)

//baliphy
//faire en sorte de le lancer plusieurs fois each x from

process baliphy{
    conda '/beegfs/data/mmorel/miniconda3/envs/baliphy'
    //publishDir "${resdir}", mode: 'copy'
    input:
    each x from {1..4}
    file sequences
    output:
    
    shell:
    '''
    bali-phy !{sequences} -S gy94[pi=f1x4] -A Codons --iter=100000000000
    '''
}
