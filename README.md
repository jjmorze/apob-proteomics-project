# apob-proteomics-project

## Plasma Proteomic Profiling of ApoB-Containing Lipoproteins in Atherosclerosis

[![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---
  
  ## Manuscript Information
  
  **Title:** Plasma Proteomic Profiling Reveals Distinct Roles of Apolipoprotein B-containing Lipoproteins in Atherosclerosis

**Authors:** Jakub Morze, Elias Björnson, Michael Y. Mi, Martin Adiels, Giorgio E. Melloni, Gonçalo Rosas da Silva, Andrzej Rynkiewicz, Nicholas A. Marston, Marta Guasch-Ferré, Michael Y. Tsai, Jerome I. Rotter, Robert E. Gerszten, Qi Sun, Frank B. Hu, Rikard Landberg, Göran Bergström, Chris J. Packard, Jan Borén, Clemens Wittenbecher

**Corresponding Authors:**  
  - Jakub Morze, MD PhD (jakub.morze@chalmers.se)  
- Clemens Wittenbecher, PhD (clemens.wittenbecher@chalmers.se)

**Affiliations:** Department of Life Sciences, Chalmers University of Technology; Institute of Medicine, University of Gothenburg; and collaborating institutions.

---
  
  ## Overview
  
  This repository contains the analysis code for generating all results, tables, and figures presented in the manuscript. The study combines observational and Mendelian randomization (MR) analyses in UK Biobank (n=35,269) and the Multi-Ethnic Study of Atherosclerosis (MESA; n=5,915) to identify plasma proteomic signatures of LDL, TRL, and Lp(a), and evaluates their association with incident coronary artery disease (CAD).

---
  
  ## Table of Contents
  
  1. [Key Findings](#key-findings)
    2. [Analytical Workflow](#analytical-workflow)
      3. [Data Availability](#data-availability)
        4. [Contact](#contact)
          5. [License](#license)
            
            ---
              
              ## Key Findings
              
              - **Proteomic signatures identified:** 30 proteins for LDL, 471 for TRL, and 53 for Lp(a)
            - **Distinct pathway enrichment:** LDL → lipid metabolism; TRL/Lp(a) → inflammatory pathways
            - **Shared proteins:** 36 proteins overlap between TRL and Lp(a) signatures (enriched in inflammation)
            - **CAD risk associations:** TRL-MPS (HR 1.16) and Lp(a)-MPS (HR 1.09) remain significant after adjustment for measured lipoprotein concentrations; LDL-MPS does not
            - **Mediation:** 62% of TRL-associated CAD risk is mediated by inflammatory protein clusters
            
            ---
              
              ## Analytical Workflow
              
              To be extended. The flowchart (overview-code-dependencies.pdf) summarizes dependencies between code notebooks, and output tables and figures. 
            
            ---
              
              ## Data Availability
              
              - **UK Biobank:** Access via [UK Biobank Access Management System](https://www.ukbiobank.ac.uk) (Applications: 570811, 53308)
            - **MESA:** Access via [MESA-NHLBI](https://www.mesa-nhlbi.org) data access procedures
            - **Summary statistics:** Provided in Supplementary Information and this repository
            
            ---
              
              ## Funding
              
              - Royal Swedish Academy of Sciences (ME2024-0009)
            - SciLifeLab & Wallenberg Data Driven Life Science Program (KAW 2020.0239)
            - Swedish Research Council (2022-01529)
            
            ---
              
              ## Contact
              
              **Jakub Morze, MD, PhD**  
              Department of Life Sciences  
            Chalmers University of Technology  
            Kemivägen 4, 41296 Gothenburg, Sweden  
            Email: jakub.morze@chalmers.se  
            GitHub: [@jjmorze](https://github.com/jjmorze)
            
            ---
              
              ## Citation
              
              ```bibtex
            @article{morze2025apob,
              title={Plasma Proteomic Profiling Reveals Distinct Roles of Apolipoprotein 
                B-containing Lipoproteins in Atherosclerosis},
              author={Morze, Jakub and Bj{\"o}rnson, Elias and Mi, Michael Y. and others},
  journal={[Journal]},
  year={2025}
}
```

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

*Code repository: https://github.com/jjmorze/apob-proteomics-project*  
*Last updated: 2025-11-27 | Version: 1.0.1*
