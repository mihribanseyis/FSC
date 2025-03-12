# **Flux-sum coupling analysis**
The Flux-sum coupling (FSC) method is based on the linear fractional programming (LFP) approach and uses the Charnes-Cooper transformation to convert it into a linear programming (LP) problem. This allows the identification of different types of metabolite coupling: directional, partial, and full coupling.
## **OS**
- Code was tested on macOS 12.7.5
## **Requirements**
- Matlab (Tested with 23.2.0.242891 (R2023b))
- [COBRA Toolbox](https://github.com/opencobra/cobratoolbox)
- Compatible genome-scale metabolic model files (e.g., iML1515, iMM904, AraCore)
## **Run Flux-sum coupling**
**To reproduce the results presented in the paper run the following scripts**
- `pFBA.m`: Calculates the flux sum of metabolites under various genetic and environmental perturbations. The metabolic model `iML1515` of *E. coli* is used in this step. Genetic and environmental perturbations are adjusted as described by Ishii et al. `met_EIs.xlsx` Includes the experimental values of metabolite concentrations and is used in this step to assess the accuracy of the flux sum in estimating these concentrations.
- `FScoupling.m`: Perform FSC analysis on metabolic models to identify coupled metabolite pairs. (iML1515, iMM904, and AraCore models are used in this step)


## REFERENCES
Ishii, Nobuyoshi, Kenji Nakahigashi, Tomoya Baba, Martin Robert, Tomoyoshi Soga, Akio Kanai, Takashi Hirasawa, et al. 2007. ‘Multiple High-Throughput Analyses Monitor the Response of E. Coli to Perturbations’. Science 316 (5824): 593–97. https://doi.org/10.1126/science.1132067.

Arnold, Anne and Zoran Nikoloski. 2014. 'Bottom-up Metabolic Reconstruction of Arabidopsis and Its Application to Determining the Metabolic Costs of Enzyme Production.' Plant Physiol 165 (3): 1380–91.

Charnes, A. and W. W. Cooper. 1962. 'Programming with Linear Fractional Functionals'. Naval Research Logistics Quarterly 9 (3‐4): 181–86.

Chung, Bevan Kai Sheng and Dong-Yup Lee. 2009. 'Flux-Sum Analysis: A Metabolite-Centric Approach for Understanding the Metabolic Network'. BMC Systems Biology 3 (1): 117.

Monk, Jonathan M. et al. 2017. 'iML1515, a Knowledgebase That Computes Escherichia Coli Traits'. Nature Biotechnology 35 (10): 904–8

