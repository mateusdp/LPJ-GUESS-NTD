LPJ-GUESS 4.1 WITH PHOSPHORUS CYCLE
========

Branch of Frankfurt-BiK-F Data Modelling center Gitlab repository

(full description to do)


PLEASE NOTE:
==============

-- 9.2022

- Phosphorus cycle balanced, including sorbed and strongly sorbed going into flux
- Balancing of labile and sorbed based on quadradic equation (extra function)
- vmax without c2 factor.
- Activating acclimated respiration gives out of memory error when running larger regions (>amazon) in Goethe HLR


-- 7.2022

- Phosphorus cycles still NOT balanced - error functions deactivated
- Balancing of labile P and sorbed P pools based on Wang et al. 2010. However using equation D8 results in inconsistencies
for soil P pools of soils type 1 and 2.
- Vmax p limiation equation based on Hidaka 2013 is Amax and conversion dividing by c2 would be necessary to convert this
into vmax, according to Haxeltine and Prentice. This is turn would make the inclusion of the tfac component redundant.
However, adding c2 to the conversion factor CN results in a very high vmax for P, and almost no P limitation globally.
Also, using tfac = 1 is what makes N limitation greater than P limitation in the Andes, which is right, but probably for
the wrong reasons.