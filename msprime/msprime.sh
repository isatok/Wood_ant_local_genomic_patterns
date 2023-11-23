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

pop_Ne_Anc=500000/2 #This seems to be the ancestor of the aqu/pol and outgroup
pop_Ne_P12_Anc=430000/2 #This is aqu/pol common ancestor; 467,232 for WSwi pol, Swi aqu
pop_Ne_H12_Anc=43/2 #NOT RELEVANT?
pop_Ne_resize_P1=310000/2 #This is aqu after aqu/pol split; 369,890 for WSwi pol, Swi aqu
pop_Ne_resize_P2=210000/2 #this is pol after aqu/pol split; 215,763 for WSwi pol, Swi aqu
pop_Ne_H1=272/2 #NOT RELEVANT?
pop_Ne_H2=168/2 #NOT RELEVANT?
pop_Ne_OG=200000/2 #OK? check?
pop_Ne_resize_P1_0=52000/2 #This is aqu after pop resize; 23,245 for WSwi pol, Swi aqu
pop_Ne_resize_P2_0=280000/2 #This is pol after pop resize; 32,522 for WSwi pol, Swi aqu
pop_Ne_P1=266/2 #WHAT'S THIS?
pop_Ne_P2=126/2 #WHAT'S THIS?
prop_P1_H=0.42
t_hyb=47 #NOT RELEVANT?
t_hyb_split=45 #NOT RELEVANT?
t_parents=225000 #This is the time of aqu/pol split; 212,802 for WSwi pol, Swi aqu
t_outgroup=2000000 # in theory, around 2e6 gens: 5Mya (Goropash. 2012) * 2.5 years / generation
t_resize=7500 #This is the timing of aqu and pol pop resize; 6,540 for WSwi pol, Swi aqu
mig_P2P1_ancestral=5.99e-6
mig_P2P1_recent=1.14e-5
mig_P1P2_recent=4.02e-6
pop_n=10
r=1e-6                 ##OK
n_blocks=100           ##OK?
block_length=1e4 # 1e4 ##OK?


###
### MODEL -----------------------------------------------------------------------------------------
###

# One parent population splits into two populations, which continue to diverge with gene flow

def sim_hybrids_shared(

pop_n, 
pop_Ne, 
pop_Ne_P1, 
pop_Ne_P2, 
pop_Ne_H1, 
pop_Ne_H2, 
pop_Ne_OG, 
pop_Ne_Anc, 
pop_Ne_P12_Anc, 
pop_Ne_H12_Anc, 
pop_Ne_resize_P1, 
pop_Ne_resize_P2, 
pop_Ne_resize_P1_0, 
pop_Ne_resize_P2_0, 
prop_P1_H, 
t_hyb, 
t_hyb_split, 
t_parents, 
t_outgroup, 
t_resize, 
mig_P2P1_ancestral, 
mig_P2P1_recent, 
mig_P1P2_recent, 
l, 
r

):
    demography = msprime.Demography()
    
    # be sure the first 5 pops are P1 P2 H1 H2 OG in this order!
    demography.add_population(name="P1", initial_size=pop_Ne_P1)
    demography.add_population(name="P2", initial_size=pop_Ne_P2)
    demography.add_population(name="H1", initial_size=pop_Ne_H1)
    demography.add_population(name="H2", initial_size=pop_Ne_H2)
    demography.add_population(name="OG", initial_size=pop_Ne_OG)
    demography.add_population(name="H12_Anc", initial_size=pop_Ne_H12_Anc)
    demography.add_population(name="Anc", initial_size=pop_Ne_Anc)
    demography.add_population(name="P12_Anc", initial_size=pop_Ne_P12_Anc)
    
    # Add recent, 2-way migration between ancestral parental populations
    demography.set_migration_rate(source="P2", dest="P1", rate=mig_P2P1_recent)
    demography.set_migration_rate(source="P1", dest="P2", rate=mig_P1P2_recent)
    
    # add split that produces two hybrid populations
    demography.add_population_split(time=t_hyb_split, derived=["H1", "H2"], ancestral="H12_Anc")
    
    # admixture event to produce hybrid ancestor population
    demography.add_admixture(time=t_hyb, derived="H12_Anc", ancestral=["P1","P2"], proportions=[prop_P1_H,1-prop_P1_H])
    
    # Reset ancestral migration from ancestral aq to ancestral pol populations
    demography.add_migration_rate_change(time=t_resize, source = "P2", dest = "P1", rate = mig_P2P1_ancestral)
    demography.add_migration_rate_change(time=t_resize, source = "P1", dest = "P2", rate = 0)
    
    # Resize parental populations
    demography.add_population_parameters_change(time=t_resize, initial_size=pop_Ne_resize_P1, population="P1")
    demography.add_population_parameters_change(time=t_resize, initial_size=pop_Ne_resize_P2, population="P2")
    
    # add population splits
    demography.add_population_split(time=t_parents, derived=["P1", "P2"], ancestral="P12_Anc")
    demography.add_population_split(time=t_outgroup, derived=["P12_Anc", "OG"], ancestral="Anc")
    
    # Resize parental populations before admixture (backwards in time)
    demography.add_population_parameters_change(time=t_hyb, initial_size=pop_Ne_resize_P1_0, population="P1")
    demography.add_population_parameters_change(time=t_hyb, initial_size=pop_Ne_resize_P2_0, population="P2")
    
    # sort events
    demography.sort_events()
    # Simulate a tree sequence
    ts = msprime.sim_ancestry(samples={"P1":pop_n, "P2":pop_n, "H1":pop_n, "H2":pop_n, "OG":1},
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
pop_Ne_H12_Anc=pop_Ne_H12_Anc, 
pop_Ne_resize_P1=pop_Ne_resize_P1, 
pop_Ne_resize_P2=pop_Ne_resize_P2, 
pop_Ne_resize_P1_0=pop_Ne_resize_P1_0, 
pop_Ne_resize_P2_0=pop_Ne_resize_P2_0, 
pop_Ne_H1=pop_Ne_H1, 
pop_Ne_H2=pop_Ne_H2, 
pop_Ne_OG=pop_Ne_OG, 
pop_Ne_P1=pop_Ne_P1, 
pop_Ne_P2=pop_Ne_P2, 
prop_P1_H=prop_P1_H, 
t_hyb=t_hyb, 
t_hyb_split=t_hyb_split, 
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

        


# Define the demographic model
demography = msprime.Demography()
demography.add_population(name="species_1", initial_size=1000)
demography.add_population(name="species_2", initial_size=1000)

# Set migration rates (gene flow) between the two species
migration_rate = 1e-4  # adjust as needed
demography.add_migration(
    source="species_1", dest="species_2", rate=migration_rate
)
demography.add_migration(
    source="species_2", dest="species_1", rate=migration_rate
)

# Set other parameters like population size changes, etc., if needed

# Simulate genetic data
tree_sequence = msprime.simulate(
    demography=demography, sample_size=20, Ne=1000, length=1e5
)

# Convert the tree sequence to VCF format
vcf_path = "simulated_data.vcf"
tree_sequence.write_vcf(vcf_path)
