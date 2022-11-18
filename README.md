--- Package under development ---

## mwss: an R package for stochastic simulation of infectious diseases spreading in healthcare systems structured as networked metapopulations

<font size="-2">
   Hammami Pachka<sup>1,2,3,*</sup>, Oodally Ajmal<sup>1,2,3,*</sup>, Reilhac Astrid<sup>4</sup>, Guérineau de Lamérie Guillaume<sup>4</sup>,  Widgren Stefan <sup>5</sup>,  Temime Laura<sup>3,6,¤</sup> and  Opatowski Lulla<sup>1,2,¤</sup></font>

<br>
<sup>1</sup> Anti-infective evasion and pharmacoepidemiology team, Université Paris-Saclay, UVSQ, Inserm, CESP,  Montigny-Le-Bretonneux, France

<sup>2</sup> Epidemiology and Modelling of Antibiotic Evasion (EMAE), Institut Pasteur, Paris, France

<sup>3</sup> Laboratoire de Modélisation,  épidémiologie et surveillance des risques sanitaires (MESuRS), Conservatoire national des arts et métiers, Paris, France

<sup>4</sup> Département d'information médicale, Centre hospitalier Guillaume Régnier, Rennes, France

<sup>5</sup> Department of Disease Control and Epidemiology, National Veterinary Institute, Uppsala, Sweden

<sup>6</sup> PACRI unit, Institut Pasteur, Conservatoire national des arts et métiers, Paris, France
<sup>7</sup>MRC Centre for Global Infectious Disease Analysis, Department of Infectious Disease Epidemiology, Imperial College London, United Kingdom

<sup>*</sup>These authors contributed equally

<sup>¤</sup>These authors contributed equally

</br>

Corresponding author: Hammami Pachka (pachka@hotmail.fr)

<!-- 
## Preprint
Preprint available at: <a href="" target="_blank"> doi: </a> 
-->
## Installation
The package can be install using 'devtools' library:

````
library(devtools)
install_github("MESuRS-Lab/mwss")
````
The companion RShiny application providing user-friendly interface to run simulations can be loaded directly from the GitHub repository (https://github.com/MESuRS-Lab/mwss-App) using 'shiny' library:
````
library(shiny) # use version >= 1.7.1
runGitHub("MESuRS-Lab/mwss-App")
````

## Main contents
This repository contains the source code for the "mwss" package developed using R-programming language.
This package allows the user to run multi-ward stochastic simulations to simulate virus transmission in the hospital setting.

````
mwss
├── R
├── data
├── man
├── DESCRIPTION
├── NAMESPACE
├── mwss.Rproj
````

- **R**
<br>  This folder contains  all the package functions

- **Data**
<br> This folder contains the toydataset to run the examples documented for each exported function.

## Further actions
- sensitivity analysis
- publications

