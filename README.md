#  Population model for the sustainable harvest of <i>Litoria raniformis</i> from wild populations
These files inform safe wild-harvest strategies for <i>Litoria raniformis</i> for conservation. These analyses underpin the manuscript "Stochastic population models guide the sustainable wild harvest of a threatened amphibian, <i>Litoria raniformis</i>"

<strong>AUTHOR</strong>: Rupert Mathwin

<strong>CONTACT</strong>: rupert.mathwin.ecology@gmail.com

<strong>URL</strong>: http://GlobalEcologyFlinders.com

<strong>INSTITUTION</strong>: Flinders University

<strong>INSTITUTION</strong>: Rupert.Mathwin.Ecology

<strong>RELEASE DATE</strong>: October 2022

R code accompanies article: 

<strong>Mathwin, R, Wassens, S, Turner, A, Heard, G, </strong> and <strong>Bradshaw, CJA</strong> Modelling effects of water regulation on the population viability of threatened apmphibians. <i>in review</i>

<strong>AIM</strong>: Stochastic population models assess the population viability of three contrasting <i>Litoria raniformis</i> populations under different harvest strategies. Harvesting strategies are: harvest as eggs, larvae or adults, and at collection proportions from 0 - 100% in 5% increments.  

Repository includes the following files:
- '<a href="https://github.com/cjabradshaw/MegafaunaSusceptibility/blob/master/matrixOperators.r">matrixOperators.R</a>' — functions to manipulate matrix models
https://github.com/RupertLovesEcology/Sustainable-Harvest-L.raniformis/blob/main/Sustainable_Harvest_Lraniformis_V10.R
- '<a href="https://github.com/RRupertLovesEcology/Sustainable-Harvest-L.raniformis/blob/main/Sustainable_Harvest_Lraniformis_V10.R">Sustainable_Harvest_Lraniformis_V10.R</a>' — R #code to run all population models and create display items.
- '<a href="https://github.com/RupertLovesEcology/RiverRegulation_Frog_PopModel/blob/main/startPopsEpsom2.csv">startPopsEpsom2.csv</a>' — 10000 burnt-in starting population for the Bendigo Water Treatment Plant.
- '<a href="https://github.com/RupertLovesEcology/RiverRegulation_Frog_PopModel/blob/main/startPopsNapNap2.csv">startPopsNapNap2.csv</a>' — 10000 burnt-in starting population for Nap Nap Waterhole.
- '<a href="https://github.com/RupertLovesEcology/RiverRegulation_Frog_PopModel/blob/main/startPopsHogwash2.csv">startPopsHogwash2.csv</a>' — 10000 burnt-in starting population for Hogwash Bend.
- '<a href="https://github.com/RupertLovesEcology/RiverRegulation_Frog_PopModel/blob/main/wetdryEpsom.csv">wetdryEpsom.csv</a>' — 10000 simulated hydrological centuries expressed as 'wet' and 'dry' sequences. wet = breeding, dry = no breeding.  For the Bendigo Water Treatment Plant models.
- '<a href="https://github.com/RupertLovesEcology/RiverRegulation_Frog_PopModel/blob/main/wetdryNapNap.csv">wetdryNapNap.csv</a>' — 10000 simulated hydrological centuries expressed as 'wet' and 'dry' sequences. wet = breeding, dry = no breeding.  For the Nap Nap Waterhole models.
- '<a href="https://github.com/RupertLovesEcology/RiverRegulation_Frog_PopModel/blob/main/wetdryHogwash.csv">wetdryHogwash.csv</a>' — 10000 simulated hydrological centuries expressed as 'wet' and 'dry' sequences. wet = breeding, dry = no breeding.  For the Hogwash Bend models.
