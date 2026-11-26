# Benthic prey contributions to Atlantic walrus diet at Southampton Island, Hudson Bay: a Bayesian stable isotope mixing model

**Author:** Bora Deniz Yılmaz

**Full report:** [walrus_mixsiar_report.pdf](walrus_mixsiar_report.pdf) · **Supplement:** [SUPPLEMENT.pdf](SUPPLEMENT.pdf)

## Abstract

Atlantic walrus (*Odobenus rosmarus rosmarus*) are large marine mammals that link seafloor production to higher trophic levels in Arctic food webs by benthic foraging. I estimated the contributions of four benthic invertebrate prey, representing different feeding guilds, to walrus diet at Southampton Island (Hudson Bay). I used published δ¹³C and δ¹⁵N values (Amiraux et al. 2023) in a Bayesian mixing model (MixSIAR). Prey were the bivalve *Hiatella arctica*, the polychaete *Nereis pelagica*, the sea cucumber *Molpadia* sp. and the gastropod *Buccinum* sp. Models were fitted with two trophic discrimination factors (Δ¹⁵N = 2.4‰ and 3.4‰). *Hiatella* was the largest single prey in both (posterior medians 53% and 68%; probability of being the largest source 0.84 and 0.96). Because *Hiatella* and *Nereis* were only partly separable (around -0.7 in pairwise comparisons), they were also combined; together they made up 72% (95% CI 48–94%) and 80% (95% CI 58–96%), depending on the TDF value used, of the diet, while sea cucumbers and whelks each contributed medians of 13% or less. Diet estimates were unaffected by how variation among individual walrus was modelled i.e whether a summary mean was used or synthetic individual data points were used. The Δ¹⁵N = 3.4‰ scenario fits previous evidence more closely, and under it *Hiatella* alone very probably formed the majority of the diet (P = 0.91). These results agree independently with stomach-content evidence that walrus rely mainly on bivalves.

## Repository contents

| File | Contents |
|---|---|
| `walrus_mixsiar_report.pdf` | Full report |
| `SUPPLEMENT.pdf` | Supplementary figures and tables |
| `analysis/Walrus_mixing.R` | Complete analysis code (MixSIAR models, mixing polygon, combined sources, probability checks) |
| `analysis/walrus_consumer.csv` | Walrus mean (Model A) |
| `analysis/walrus_consumer_n4_synthetic.csv` | Four synthetic walrus reproducing the reported mean and SD (Model B) |
| `analysis/walrus_sources.csv` | Prey isotope means, SDs and n |
| `analysis/walrus_discrimination_main.csv` | TDF-H: Hobson et al. (1996) |
| `analysis/walrus_discrimination_post.csv` | TDF-P: Δ¹⁵N from Post (2002) |

## Requirements

R 4.6.2, JAGS 4.3.2, MixSIAR 3.1.12, ggplot2.

## Data source

All isotope values are from Amiraux et al. (2023), published open access under a CC BY 4.0 licence:

Amiraux R, Mundy CJ, Pierrejean M, et al. (2023) Tracing carbon flow and trophic structure of a coastal Arctic marine food web using highly branched isoprenoids and carbon, nitrogen and sulfur stable isotopes. *Ecological Indicators* 147: 109938. https://doi.org/10.1016/j.ecolind.2023.109938
