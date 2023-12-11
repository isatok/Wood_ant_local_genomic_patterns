#### This script is for running demography_debug() for the simulated demographies to make sure they are correct.

## The model parameter values are from Portinha et. al 2022 (https://doi.org/10.1111/mec.16481), Table S5 & S6, "Sympatry (asymmetrical mig.)".

conda activate msprime-env
python
import msprime
import numpy as np


### For Finnish F. aquilonia & Finnish F. polyctena ### --------------------------------------------------------------------------------------

demography = msprime.Demography()
pop_Ne_Anc=500000/2            #outgroup & aqu/pol common ancestor Ne
pop_Ne_OG=200000/2             #outgroup Ne
pop_Ne_P12_Anc=430000/2        #aqu & pol common ancestor Ne (backw in time: before outgroup merge)
pop_Ne_resize_P1=310000/2      #aqu Ne, after resize/before clade merge (backw in time)
pop_Ne_resize_P2=210000/2      #pol Ne, after resize/before clade merge (backw in time)
pop_Ne_P1=52000/2              #the initial aqu Ne
pop_Ne_P2=280000/2             #the initial pol Ne
t_parents=225000               #timing of aqu & pol split (i.e. merge backw in time)
t_outgroup=2000000             #timing of outgroup & aqu/pol ancestor split (i.e. merge backw in time). In theory, around 2e6 gens: 5Mya (Goropash. 2012) / 2.5 years per generation
t_resize=7500                  #timing of aqu & pol resize
mig_P2P1_ancestral=5.99e-6     #ancestral migration rate from aqu to pol (in msprime's backwards-world the migrating lineages go from pol to aqu)
mig_P2P1_recent=1.14e-5        #recent migration rate from aqu to pol (in msprime's backwards-world the migrating lineages go from pol to aqu)
mig_P1P2_recent=4.02e-6        #recent migration rate from pol to aqu (in msprime's backwards-world the migrating lineages go from aqu to pol)

demography.add_population(name="P1", initial_size=pop_Ne_P1)
demography.add_population(name="P2", initial_size=pop_Ne_P2)
demography.add_population(name="OG", initial_size=pop_Ne_OG)
demography.add_population(name="Anc", initial_size=pop_Ne_Anc)
demography.add_population(name="P12_Anc", initial_size=pop_Ne_P12_Anc)
demography.set_migration_rate(source="P2", dest="P1", rate=mig_P2P1_recent)
demography.set_migration_rate(source="P1", dest="P2", rate=mig_P1P2_recent)
demography.add_migration_rate_change(time=t_resize, source = "P2", dest = "P1", rate = mig_P2P1_ancestral)
demography.add_migration_rate_change(time=t_resize, source = "P1", dest = "P2", rate = 0)
demography.add_population_parameters_change(time=t_resize, initial_size=pop_Ne_resize_P1, population="P1")
demography.add_population_parameters_change(time=t_resize, initial_size=pop_Ne_resize_P2, population="P2")
demography.add_population_split(time=t_parents, derived=["P1", "P2"], ancestral="P12_Anc")
demography.add_population_split(time=t_outgroup, derived=["P12_Anc", "OG"], ancestral="Anc")
demography.sort_events()

demography.debug()
print(demography.debug())

DemographyDebugger
╠════════════════════════════════════╗
║ Epoch[0]: [0, 7.5e+03) generations ║
╠════════════════════════════════════╝
╟    Populations (total=5 active=3)
║    ┌─────────────────────────────────────────────────────────────────────┐
║    │    │      start│        end│growth_rate  │    P1    │    P2    │ OG │
║    ├─────────────────────────────────────────────────────────────────────┤
║    │  P1│    26000.0│    26000.0│ 0           │    0     │ 4.02e-06 │ 0  │
║    │  P2│   140000.0│   140000.0│ 0           │ 1.14e-05 │    0     │ 0  │
║    │  OG│   100000.0│   100000.0│ 0           │    0     │    0     │ 0  │
║    └─────────────────────────────────────────────────────────────────────┘
╟    Events @ generation 7.5e+03
║    ┌────────────────────────────────────────────────────────────────────────────────────────┐
║    │  time│type            │parameters             │effect                                  │
║    ├────────────────────────────────────────────────────────────────────────────────────────┤
║    │  7500│Migration rate  │source=P2, dest=P1,    │Backwards-time migration rate from P2   │
║    │      │change          │rate=5.99e-06          │to P1 → 5.99e-06                        │
║    │┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈│
║    │  7500│Migration rate  │source=P1, dest=P2,    │Backwards-time migration rate from P1   │
║    │      │change          │rate=0                 │to P2 → 0                               │
║    │┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈│
║    │  7500│Population      │population=P1,         │initial_size → 1.6e+05 for population   │
║    │      │parameter       │initial_size=155000.0  │P1                                      │
║    │      │change          │                       │                                        │
║    │┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈│
║    │  7500│Population      │population=P2,         │initial_size → 1e+05 for population P2  │
║    │      │parameter       │initial_size=105000.0  │                                        │
║    │      │change          │                       │                                        │
║    └────────────────────────────────────────────────────────────────────────────────────────┘
╠═══════════════════════════════════════════╗
║ Epoch[1]: [7.5e+03, 2.25e+05) generations ║
╠═══════════════════════════════════════════╝
╟    Populations (total=5 active=3)
║    ┌───────────────────────────────────────────────────────────────┐
║    │    │      start│        end│growth_rate  │    P1    │ P2 │ OG │
║    ├───────────────────────────────────────────────────────────────┤
║    │  P1│   155000.0│   155000.0│ 0           │    0     │ 0  │ 0  │
║    │  P2│   105000.0│   105000.0│ 0           │ 5.99e-06 │ 0  │ 0  │
║    │  OG│   100000.0│   100000.0│ 0           │    0     │ 0  │ 0  │
║    └───────────────────────────────────────────────────────────────┘
╟    Events @ generation 2.25e+05
║    ┌──────────────────────────────────────────────────────────────────────────────────┐
║    │      time│type        │parameters         │effect                                │
║    ├──────────────────────────────────────────────────────────────────────────────────┤
║    │  2.25e+05│Population  │derived=[P1, P2],  │Moves all lineages from derived       │
║    │          │Split       │ancestral=P12_Anc  │populations 'P1' and 'P2' to the      │
║    │          │            │                   │ancestral 'P12_Anc' population. Also  │
║    │          │            │                   │set the derived populations to        │
║    │          │            │                   │inactive, and all migration rates to  │
║    │          │            │                   │and from the derived populations to   │
║    │          │            │                   │zero.                                 │
║    └──────────────────────────────────────────────────────────────────────────────────┘
╠═════════════════════════════════════════╗
║ Epoch[2]: [2.25e+05, 2e+06) generations ║
╠═════════════════════════════════════════╝
╟    Populations (total=5 active=2)
║    ┌──────────────────────────────────────────────────────────────┐
║    │         │      start│        end│growth_rate  │ OG │ P12_Anc │
║    ├──────────────────────────────────────────────────────────────┤
║    │       OG│   100000.0│   100000.0│ 0           │ 0  │    0    │
║    │  P12_Anc│   215000.0│   215000.0│ 0           │ 0  │    0    │
║    └──────────────────────────────────────────────────────────────┘
╟    Events @ generation 2e+06
║    ┌─────────────────────────────────────────────────────────────────────────────────────┐
║    │   time│type        │parameters              │effect                                 │
║    ├─────────────────────────────────────────────────────────────────────────────────────┤
║    │  2e+06│Population  │derived=[P12_Anc, OG],  │Moves all lineages from derived        │
║    │       │Split       │ancestral=Anc           │populations 'P12_Anc' and 'OG' to the  │
║    │       │            │                        │ancestral 'Anc' population. Also set   │
║    │       │            │                        │the derived populations to inactive,   │
║    │       │            │                        │and all migration rates to and from    │
║    │       │            │                        │the derived populations to zero.       │
║    └─────────────────────────────────────────────────────────────────────────────────────┘
╠════════════════════════════════════╗
║ Epoch[3]: [2e+06, inf) generations ║
╠════════════════════════════════════╝
╟    Populations (total=5 active=1)
║    ┌───────────────────────────────────────────┐
║    │     │      start│        end│growth_rate  │
║    ├───────────────────────────────────────────┤
║    │  Anc│   250000.0│   250000.0│ 0           │
║    └───────────────────────────────────────────┘


### For Swiss F. aquilonia & West Swiss F. polyctena ### --------------------------------------------------------------------------------------

demography = msprime.Demography()

pop_Ne_Anc=500000/2            #outgroup & aqu/pol common ancestor Ne
pop_Ne_OG=200000/2             #outgroup Ne
pop_Ne_P12_Anc=467000/2        #aqu & pol common ancestor Ne (backw in time: before outgroup merge)
pop_Ne_resize_P1=370000/2      #aqu Ne, after resize/before clade merge (backw in time)
pop_Ne_resize_P2=216000/2      #pol Ne, after resize/before clade merge (backw in time)
pop_Ne_P1=23000/2              #the initial aqu Ne
pop_Ne_P2=33000/2              #the initial pol Ne
t_parents=213000               #timing of aqu & pol split (i.e. merge backw in time)
t_outgroup=2000000             #timing of outgroup & aqu/pol ancestor split (i.e. merge backw in time). In theory, around 2e6 gens: 5Mya (Goropash. 2012) / 2.5 years per generation
t_resize=6500                  #timing of aqu & pol resize
mig_P2P1_ancestral=2.36e-6     #ancestral migration rate from aqu to pol (in msprime's backwards-world the migrating lineages go from pol to aqu).
mig_P2P1_recent=3.81e-6        #recent migration rate from aqu to pol (in msprime's backwards-world the migrating lineages go from pol to aqu).
mig_P1P2_recent=3.35e-12       #recent migration rate from pol to aqu (in msprime's backwards-world the migrating lineages go from aqu to pol).


demography.add_population(name="P1", initial_size=pop_Ne_P1)
demography.add_population(name="P2", initial_size=pop_Ne_P2)
demography.add_population(name="OG", initial_size=pop_Ne_OG)
demography.add_population(name="Anc", initial_size=pop_Ne_Anc)
demography.add_population(name="P12_Anc", initial_size=pop_Ne_P12_Anc)
demography.set_migration_rate(source="P2", dest="P1", rate=mig_P2P1_recent)
demography.set_migration_rate(source="P1", dest="P2", rate=mig_P1P2_recent)
demography.add_migration_rate_change(time=t_resize, source = "P2", dest = "P1", rate = mig_P2P1_ancestral)
demography.add_migration_rate_change(time=t_resize, source = "P1", dest = "P2", rate = 0)
demography.add_population_parameters_change(time=t_resize, initial_size=pop_Ne_resize_P1, population="P1")
demography.add_population_parameters_change(time=t_resize, initial_size=pop_Ne_resize_P2, population="P2")
demography.add_population_split(time=t_parents, derived=["P1", "P2"], ancestral="P12_Anc")
demography.add_population_split(time=t_outgroup, derived=["P12_Anc", "OG"], ancestral="Anc")
demography.sort_events()

demography.debug()
print(demography.debug())
╠════════════════════════════════════╗
║ Epoch[0]: [0, 6.5e+03) generations ║
╠════════════════════════════════════╝
╟    Populations (total=5 active=3)
║    ┌─────────────────────────────────────────────────────────────────────┐
║    │    │      start│        end│growth_rate  │    P1    │    P2    │ OG │
║    ├─────────────────────────────────────────────────────────────────────┤
║    │  P1│    11500.0│    11500.0│ 0           │    0     │ 3.35e-12 │ 0  │
║    │  P2│    16500.0│    16500.0│ 0           │ 3.81e-06 │    0     │ 0  │
║    │  OG│   100000.0│   100000.0│ 0           │    0     │    0     │ 0  │
║    └─────────────────────────────────────────────────────────────────────┘
╟    Events @ generation 6.5e+03
║    ┌───────────────────────────────────────────────────────────────────────────────────────┐
║    │  time│type            │parameters             │effect                                 │
║    ├───────────────────────────────────────────────────────────────────────────────────────┤
║    │  6500│Migration rate  │source=P2, dest=P1,    │Backwards-time migration rate from P2  │
║    │      │change          │rate=2.36e-06          │to P1 → 2.36e-06                       │
║    │┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈│
║    │  6500│Migration rate  │source=P1, dest=P2,    │Backwards-time migration rate from P1  │
║    │      │change          │rate=0                 │to P2 → 0                              │
║    │┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈│
║    │  6500│Population      │population=P1,         │initial_size → 1.8e+05 for population  │
║    │      │parameter       │initial_size=185000.0  │P1                                     │
║    │      │change          │                       │                                       │
║    │┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈│
║    │  6500│Population      │population=P2,         │initial_size → 1.1e+05 for population  │
║    │      │parameter       │initial_size=108000.0  │P2                                     │
║    │      │change          │                       │                                       │
║    └───────────────────────────────────────────────────────────────────────────────────────┘
╠═══════════════════════════════════════════╗
║ Epoch[1]: [6.5e+03, 2.13e+05) generations ║
╠═══════════════════════════════════════════╝
╟    Populations (total=5 active=3)
║    ┌───────────────────────────────────────────────────────────────┐
║    │    │      start│        end│growth_rate  │    P1    │ P2 │ OG │
║    ├───────────────────────────────────────────────────────────────┤
║    │  P1│   185000.0│   185000.0│ 0           │    0     │ 0  │ 0  │
║    │  P2│   108000.0│   108000.0│ 0           │ 2.36e-06 │ 0  │ 0  │
║    │  OG│   100000.0│   100000.0│ 0           │    0     │ 0  │ 0  │
║    └───────────────────────────────────────────────────────────────┘
╟    Events @ generation 2.13e+05
║    ┌──────────────────────────────────────────────────────────────────────────────────┐
║    │      time│type        │parameters         │effect                                │
║    ├──────────────────────────────────────────────────────────────────────────────────┤
║    │  2.13e+05│Population  │derived=[P1, P2],  │Moves all lineages from derived       │
║    │          │Split       │ancestral=P12_Anc  │populations 'P1' and 'P2' to the      │
║    │          │            │                   │ancestral 'P12_Anc' population. Also  │
║    │          │            │                   │set the derived populations to        │
║    │          │            │                   │inactive, and all migration rates to  │
║    │          │            │                   │and from the derived populations to   │
║    │          │            │                   │zero.                                 │
║    └──────────────────────────────────────────────────────────────────────────────────┘
╠═════════════════════════════════════════╗
║ Epoch[2]: [2.13e+05, 2e+06) generations ║
╠═════════════════════════════════════════╝
╟    Populations (total=5 active=2)
║    ┌──────────────────────────────────────────────────────────────┐
║    │         │      start│        end│growth_rate  │ OG │ P12_Anc │
║    ├──────────────────────────────────────────────────────────────┤
║    │       OG│   100000.0│   100000.0│ 0           │ 0  │    0    │
║    │  P12_Anc│   233500.0│   233500.0│ 0           │ 0  │    0    │
║    └──────────────────────────────────────────────────────────────┘
╟    Events @ generation 2e+06
║    ┌─────────────────────────────────────────────────────────────────────────────────────┐
║    │   time│type        │parameters              │effect                                 │
║    ├─────────────────────────────────────────────────────────────────────────────────────┤
║    │  2e+06│Population  │derived=[P12_Anc, OG],  │Moves all lineages from derived        │
║    │       │Split       │ancestral=Anc           │populations 'P12_Anc' and 'OG' to the  │
║    │       │            │                        │ancestral 'Anc' population. Also set   │
║    │       │            │                        │the derived populations to inactive,   │
║    │       │            │                        │and all migration rates to and from    │
║    │       │            │                        │the derived populations to zero.       │
║    └─────────────────────────────────────────────────────────────────────────────────────┘
╠════════════════════════════════════╗
║ Epoch[3]: [2e+06, inf) generations ║
╠════════════════════════════════════╝
╟    Populations (total=5 active=1)
║    ┌───────────────────────────────────────────┐
║    │     │      start│        end│growth_rate  │
║    ├───────────────────────────────────────────┤
║    │  Anc│   250000.0│   250000.0│ 0           │
║    └───────────────────────────────────────────┘

