import numpy as np
import msprime
import gzip
import sys


#### Model #### ---------------------------------------------------------------------------------


# An ancestral population splits into two populations ("F. aquilonia" (P1 "aqu"), and "F. polyctena" (P2 "pol")), which continue to diverge with gene flow

def aqu_pol(
    pop_n,                 #number of analysed diploid samples per group (in this case each group is a species)  #####<-----IS THIS CORRECT?
    pop_Ne_OG,             #outgroup Ne
    pop_Ne_Anc,            #outgroup & aqu/pol common ancestor Ne
    pop_Ne_P12_Anc,        #aqu & pol common ancestor Ne
    pop_Ne_resize_P1,      #aqu Ne after resize (backw in time)
    pop_Ne_resize_P2,      #pol Ne after resize (backw in time)
    pop_Ne_P1,             #the initial aqu Ne
    pop_Ne_P2,             #the initial pol Ne
    t_parents,             #timing of aqu & pol split
    t_outgroup,            #timing of outgroup & aqu/pol split
    t_resize,              #timing of aqu & pol resize
    mig_P2P1_ancestral,    #ancestral migration rate from aqu to pol (in msprime's backwards-world the migrating lineages go from pol to aqu)
    mig_P2P1_recent,       #recent migration rate from aqu to pol (in msprime's backwards-world the migrating lineages go from pol to aqu)
    mig_P1P2_recent,       #recent migration rate from pol to aqu (in msprime's backwards-world the migrating lineages go from aqu to pol)
    l,                     #length of genomic blocks for which window-based stats will be calculated (each block = one window) #####<-----IS THIS CORRECT?
    r                      #recombination rate
):
    demography = msprime.Demography()
    
    # be sure the first 5 pops are P1 P2 H1 H2 OG in this order! #####<-----What is the meaning of this? Can I take H pops off from the middle without changing things elsewhere?
    
    demography.add_population(name="P1", initial_size=pop_Ne_P1)
    demography.add_population(name="P2", initial_size=pop_Ne_P2)
    demography.add_population(name="OG", initial_size=pop_Ne_OG)
    demography.add_population(name="Anc", initial_size=pop_Ne_Anc)
    demography.add_population(name="P12_Anc", initial_size=pop_Ne_P12_Anc)
    
    # Add recent, 2-way migration between ancestral parental populations
    demography.set_migration_rate(source="P2", dest="P1", rate=mig_P2P1_recent)
    demography.set_migration_rate(source="P1", dest="P2", rate=mig_P1P2_recent)
    
    # Reset ancestral migration from ancestral aq to ancestral pol populations
    demography.add_migration_rate_change(time=t_resize, source = "P2", dest = "P1", rate = mig_P2P1_ancestral)
    demography.add_migration_rate_change(time=t_resize, source = "P1", dest = "P2", rate = 0)
    
    # Resize parental populations
    demography.add_population_parameters_change(time=t_resize, initial_size=pop_Ne_resize_P1, population="P1")
    demography.add_population_parameters_change(time=t_resize, initial_size=pop_Ne_resize_P2, population="P2")
    
    # add population splits
    demography.add_population_split(time=t_parents, derived=["P1", "P2"], ancestral="P12_Anc")
    demography.add_population_split(time=t_outgroup, derived=["P12_Anc", "OG"], ancestral="Anc")
    
    # sort events
    demography.sort_events()
    # Simulate a tree sequence
    ts = msprime.sim_ancestry(samples={"P1":pop_n, "P2":pop_n, "OG":1},
                              demography=demography, ploidy = 2, sequence_length = l, recombination_rate = r)
    return(ts)


#### Parameters #### ---------------------------------------------------------------------------------


# An ancestral population splits into two populations ("F. aquilonia" (P1 "aqu"), and "F. polyctena" (P2 "pol")), which continue to diverge with gene flow

pop_Ne_Anc=500000/2            #outgroup & aqu/pol common ancestor Ne
pop_Ne_OG=200000/2             #outgroup Ne
pop_Ne_P12_Anc=430000/2        #aqu & pol common ancestor Ne (backw in time: before outgroup merge); 467,232 for WSwi pol, Swi aqu

pop_Ne_resize_P1=310000/2      #aqu Ne, after resize/before clade merge (backw in time); 369,890 for WSwi pol, Swi aqu
pop_Ne_resize_P2=210000/2      #pol Ne, after resize/before clade merge (backw in time); 215,763 for WSwi pol, Swi aqu

pop_Ne_P1=52000/2              #the initial aqu Ne; 23,245 for WSwi pol, Swi aqu
pop_Ne_P2=280000/2             #the initial pol Ne; 32,522 for WSwi pol, Swi aqu

t_parents=225000               #timing of aqu & pol split (i.e. merge backw in time); 212,802 for WSwi pol, Swi aqu
t_outgroup=2000000             #timing of outgroup & aqu/pol ancestor split (i.e. merge backw in time). In theory, around 2e6 gens: 5Mya (Goropash. 2012) / 2.5 years per generation

t_resize=7500                  #timing of aqu & pol resize; 6,540 for WSwi pol, Swi aqu

mig_P2P1_ancestral=5.99e-6     #ancestral migration rate from aqu to pol (in msprime's backwards-world the migrating lineages go from pol to aqu)
mig_P2P1_recent=1.14e-5        #recent migration rate from aqu to pol (in msprime's backwards-world the migrating lineages go from pol to aqu)
mig_P1P2_recent=4.02e-6        #recent migration rate from pol to aqu (in msprime's backwards-world the migrating lineages go from aqu to pol)

pop_n=10                       #number of analysed diploid samples per group (in this case each group is a species)
r=1e-6                         #recombination rate in centimorgans per basepair (equals to 1cM/Mb)
n_blocks=100                   #number of genomic blocks for which window-based stats will be calculated (each block = one window)
block_length=1e4               #length of genomic blocks for which window-based stats will be calculated (each block = one window)



#### Simulation #### ---------------------------------------------------------------------------------


# run simulations to create tree sequence objects

ts_blocks = [aqu_pol(

pop_n,  
pop_Ne_Anc=pop_Ne_Anc, 
pop_Ne_P12_Anc=pop_Ne_P12_Anc,  
pop_Ne_resize_P1=pop_Ne_resize_P1, 
pop_Ne_resize_P2=pop_Ne_resize_P2, 
pop_Ne_P1=pop_Ne_P1, 
pop_Ne_P2=pop_Ne_P2, 
pop_Ne_OG=pop_Ne_OG, 
t_parents=t_parents, 
t_outgroup=t_outgroup, 
t_resize=t_resize, 
mig_P2P1_ancestral=mig_P2P1_ancestral, 
mig_P2P1_recent=mig_P2P1_recent, 
mig_P1P2_recent=mig_P1P2_recent, 
l=block_length, r=r) 

for i in range(n_blocks)]




#### VCF #### ---------------------------------------------------------------------------------


# add mutations

ts_blocks_mutated = [None]*n_blocks
for i in range(n_blocks):
    ts_blocks_mutated[i] = msprime.sim_mutations(ts_blocks[i], rate=3.5e-9)
    
# write VCF file ##########MIKÄ WT?????#########

with gzip.open("../output.vcf.gz", "wt") as vcf_file:
    for i in range(n_blocks):
        ts_blocks_mutated[i].write_vcf(vcf_file)



