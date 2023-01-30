LPJ-GUESS 4.1 WITH PHOSPHORUS CYCLE, TRAIT VARIATION AND MYCORRHIZAE (CF input implemented)
========

Branch of Frankfurt-BiK-F Data Modelling center Gitlab repository

(full description to do)


PLEASE NOTE:
==============

-- 30.1.2023

To add, turnover tradeoff with a maximum of 1 (limitation of yearly growth/dynamics)


-- 24.1.2023

Added differential NO3 or NH4 uptake for roots and mycorrhiza


-- 21.10.2022

Added sla variation - changes results but no large difference
Added cton and ctop tradeoffs - reduces biomass and especially cton and ctop ratios.
Added wsg variation - increases biomass in relation to previous implementation, but still relatively low.
Added output variables for SLA and WSG
obs. cton and ctop tradeoffs for Ecuador when added (indiv.ctonp_tradeoffs()) cause lower cton and ctop than 
standard and a significant reduction in cmass, lai and npp.
Using the correct kmp values improves biomass. Maybe the model needs P myco uptake.


-- 20.10.2022

Added variables to individual, maybe issues with BLAZE, transfer litter, canexch.cpp line 2931 because still with pft cton variables.


-- 18.10.2022

Started development from phosphorous branch.