//path_file = "/beegfs/data/mmorel/catarina_gls/"
path_file = "/beegfs/data/mmorel/Infer_Anth/align/"


params.sequences = path_file + "seed_plants_cycas_domains_aa.fa"
//"subset_seed_plants.aa.fa"
//params.resdir= path_file + "results_baliphy_complete/" //do not forget to change the resdir

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
    bali-phy !{sequences} -S 'lg08+Rates.gamma+inv' --iter=100000000000
    '''
}
