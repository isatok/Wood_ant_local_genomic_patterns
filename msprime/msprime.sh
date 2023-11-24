test script

import numpy as np
import msprime
import twisst
import gzip
import sys

###
### PARAMETERS ------------------------------------------------------------------------------------
###

# One parent population splits into two populations, which continue to diverge with gene flow

pop_Ne_Anc=500000/2 #Ne for the common ancestor of aqu, pol, and the outgroup
pop_Ne_OG=200000/2 #Ne for the outgroup
pop_Ne_P12_Anc=430000/2 #Ne for aqu and pol common ancestor (backw in time: before outgroup merge); 467,232 for WSwi pol, Swi aqu

pop_Ne_resize_P1=310000/2 #New Ne for aqu when the population was resized (backw in time: before clade merge); 369,890 for WSwi pol, Swi aqu
pop_Ne_resize_P2=210000/2 #New Ne for pol when the population was resized (backw in time: before clade merge); 215,763 for WSwi pol, Swi aqu

pop_Ne_resize_P1_0=52000/2 #Ne for aqu from present to first resize (backw in time); 23,245 for WSwi pol, Swi aqu
pop_Ne_resize_P2_0=280000/2 #Ne for pol from present to first resize (backw in time); 32,522 for WSwi pol, Swi aqu

t_parents=225000 #Timing of aqu/pol lineage merge; 212,802 for WSwi pol, Swi aqu
t_outgroup=2000000 #Timing of outgroup lineage merge to aqu/pol ancestor. In theory, around 2e6 gens: 5Mya (Goropash. 2012) / 2.5 years per generation

t_resize=7500 #Timing of aqu and pol pop resize (after this, backw in time, e.g. Ne pol = pop_Ne_resize_P2) ; 6,540 for WSwi pol, Swi aqu

mig_P2P1_ancestral=5.99e-6 #0.000 005 99 Migration from aqu to pol after resize (i.e. before aqu & pol merge backw in time)
mig_P2P1_recent=1.14e-5    #0.000 011 4 Migration from aqu to pol recently (i.e. before aqu & pol pop resize back in time)
mig_P1P2_recent=4.02e-6 #Migration from pol to aqu recently (i.e. before aqu & pol pop resize back in time)

pop_n=10 #Number of analysed "populations"
r=1e-6 #Recombination rate
n_blocks=100 #Number of "blocks" that do not recombine with each other, within every replicate
block_length=1e4 #Length of genomic "blocks"; recombination within each 10 000 bp block is allowed

###
### MODEL -----------------------------------------------------------------------------------------
###

# One parent population splits into two populations, which continue to diverge with gene flow

def sim_hybrids_shared(pop_n, pop_Ne, pop_Ne_OG, pop_Ne_Anc, pop_Ne_P12_Anc, pop_Ne_resize_P1, pop_Ne_resize_P2, pop_Ne_resize_P1_0, 
pop_Ne_resize_P2_0, t_parents, t_outgroup, t_resize, mig_P2P1_ancestral, mig_P2P1_recent, mig_P1P2_recent, l, r):
    demography = msprime.Demography()
    
    # be sure the first 5 pops are P1 P2 H1 H2 OG in this order! ##What is the meaning of this? Can I take pops off from the middle; 
    # does it affect the order of the other pops?
    demography.add_population(name="P1", initial_size=pop_Ne_resize_P1_0)
    demography.add_population(name="P2", initial_size=pop_Ne_resize_P2_0)
    demography.add_population(name="OG", initial_size=pop_Ne_OG)
    demography.add_population(name="H12_Anc", initial_size=pop_Ne_H12_Anc)
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


###
### SIMULATE --------------------------------------------------------------------------------------
###

# run simulations to produce tree sequence objects

ts_shared_blocks = [sim_hybrids_shared(

pop_n, 
pop_Ne=1e4, 
pop_Ne_Anc=pop_Ne_Anc, 
pop_Ne_P12_Anc=pop_Ne_P12_Anc,  
pop_Ne_resize_P1=pop_Ne_resize_P1, 
pop_Ne_resize_P2=pop_Ne_resize_P2, 
pop_Ne_resize_P1_0=pop_Ne_resize_P1_0, 
pop_Ne_resize_P2_0=pop_Ne_resize_P2_0, 
pop_Ne_OG=pop_Ne_OG, 
t_parents=t_parents, 
t_outgroup=t_outgroup, 
t_resize=t_resize, 
mig_P2P1_ancestral=mig_P2P1_ancestral, 
mig_P2P1_recent=mig_P2P1_recent, 
mig_P1P2_recent=mig_P1P2_recent, 
l=block_length, r=r) 

for i in range(n_blocks)]



###
### VCF -------------------------------------------------------------------------------------------
###

# add mutations

ts_shared_blocks_mutated = [None]*n_blocks
for i in range(n_blocks):
    ts_shared_blocks_mutated[i] = msprime.sim_mutations(ts_shared_blocks[i], rate=3.5e-9) ##OK
    
# write VCF file

with gzip.open("../output.vcf.gz", "wt") as vcf_file:
    for i in range(n_blocks):
        ts_shared_blocks_mutated[i].write_vcf(vcf_file)



