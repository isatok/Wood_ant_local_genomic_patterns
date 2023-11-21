test script

import msprime

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
