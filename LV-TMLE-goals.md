# Goals of the Latent-Variable Targeted Maximum Likelihood Estimation (LV-TMLE)

### This document outlines to broad goals of this research project.

### Overarching Goals:
1. To integrate latent variables from the structural equation modeling tradition with targeted maximum likelihood estimation (LV-TMLE).
2. To successfully conduct a simulation study that compares LV-TMLE with suitable alternatives.
3. To demonstrate that LV-TMLE is a viable alternative in realistic conditions for psychologists to use. For LV-TMLE to be viable, it should be better than the alternatives we compare it to in our simulation. This means that it should perform better in most conditions but LV-TMLE does not have to perform better in all conditions.

### Constraints
- Manuscript writing should be in a style and format that is suitable to publish in the journal, "Multivariate Behavioral Research".
- Replications for each condition within the simulation should be conducted 1000 times.
- The full simulation should be run on the Bouchet Yale Computing Cluster, given the computational resources required.
- When writing the manuscript, it is unacceptable to have hallucinated references. When writing and citing authors, extra care must be taken to make sure the reference exists. For example, you should use the "deep research" feature and any other tools at your disposal to ensure all references and citations are linked to real papers.
- Statistical rigor and correct derivations are essential. Break down any proofs into small sub-tasks to avoid any mistakes, especially when derivations require many steps. 
- Consult the "quant-writing-skills.md" file before doing any manuscript writing.

### Important Sub-goals
1. LV-TMLE should have better coverage rates in most simulation conditions than the alternatives, especially when model misspecification exists.
2. To integrate latent variables with TMLE, multiple methodological innovations are required. For example, using variants of TMLE that can accommodate continuous exposure variables, techniques to draw better samples from the posterior distribution etc.). These innovations and features should be documented so they can possibly be highlighted in the eventual manuscript.
3. To have reproducible and well-documented code that human researchers can read and easily understand.
4. The bias and RMSE should be lower in the LV-TMLE than the other models. While LV-TMLE cannot outperform alternatives such as SEM in all conditions, LV-TMLE should do better in conditions that are realistic (e.g., realistic sample sizes in psychology) such that applied researchers are actually incentivized to use LV-TMLE.
5. The immediate goal is not create an R package. However, any code you write should keep in mind that I eventually want to make LV-TMLE available for applied researchers to use. 

