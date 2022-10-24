## Research compendium for 'A 2D geometric morphometric assessment of chrono-cultural trends in osseous barbed points of the European Final Palaeolithic and Early Mesolithic' 

### Compendium DOI:

[![DOI](https://zenodo.org/badge/DOI/.svg)](https://doi.org/)

The files at the URL above will generate the results as found in the publication. The files hosted at <https://github.com/yesdavid/Osseous_Barbed_Points_2D_GMM> are the development versions and may have changed since the paper was published.

### Maintainer of this repository:

[![ORCiD](https://img.shields.io/badge/ORCiD-0000--0001--7349--5401-green.svg)](http://orcid.org/0000-0001-7349-5401) David N. Matzig (<david.matzig@cas.au.dk>) 

### Published in:

Tsirintoulaki, K., Matzig, D.N., Riede, F. (submitted). A 2D geometric morphometric assessment of chrono-cultural trends in osseous barbed points of the European Final Palaeolithic and Early Mesolithic. _Open Archaeology_. [![DOI](https://zenodo.org/badge/DOI/.svg)](https://doi.org/).

### Abstract:

Studies on prehistoric osseous barbed points have relied heavily on typology in linking presumed types to broader techno-complexes and for chronological inference. The accretion of both new finds and of radiocarbon dates obtained directly on such artefacts, however, has revealed that i) shape variability defies neat typological divisions, and that ii) chronological inferences based on typology often fail. To further query these issues and to better understand the design choices and cultural evolutionary dynamics within this artefact class, we present a 2D whole-outline geometric morphometric analysis of 50 directly dated Late Pleistocene and Early Holocene osseous barbed points from primarily northern and western Europe. The results indicate that a) different components (tip, base, barbs) of these artefacts were subject to varying design constraints and that b) there is no clear-cut distinction between Final Palaeolithic and Mesolithic point traditions. Different techno-functional components evolved at various rates while specimens assigned to the same type and/or techno-complex only occasionally are morphologically similar. The results reflect a relatively low level of normativity for this artefact class and repeated convergence on similar design elements. We propose that interpretations linked to cultural dynamics, individual craft agency and repeated convergence on locally optimal designs may offer more satisfying avenues for thinking about the barbed points of this period.

### Keywords: 

Barbed points, Late Palaeolithic and Mesolithic, geometric morphometrics, cultural transmission

### Overview of contents and how to reproduce:

This repository contains code (`1_script`) data (`2_data`) and for the paper. After downloading, the results can be reproduced using ` Tsirintoulaki_et_al_Barbed_Points.Rproj` and the existing folder structure. All analyses and visualisations presented in this paper were prepared in R version 4.2.1 (2022-06-23) under Ubuntu 18.04.5 LTS (64-bit).

### Required dependencies:

As the data and code in this repository are complete and self-contained, it can be reproduced with only an R environment (tested for R v4.2.1). The necessary package dependencies are documented in the DESCRIPTION file and can be installed manually or automatically with

```
if(!require("remotes")) install.packages("remotes")
remotes::install_github("yesdavid/Osseous_Barbed_Points_2D_GMM", 
repos = "https://mran.microsoft.com/snapshot/2022-09-30")
```

This will install the relevant package dependency versions from September 2022, thanks to Microsoft's [CRAN Time Machine](https://mran.microsoft.com/timemachine).

### Licenses:

Code: MIT <http://opensource.org/licenses/MIT> year: 2022, copyright holder: David Nicolas Matzig
