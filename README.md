# Masters_Thesis-project

# Optimizing Single-Cell RNA Sequencing Integration: A Dual-Pipeline Approach
Comparative Benchmarking of Biological Conservation, Batch Correction, and Computational Efficiency

## Introduction
This repository contains the code and data used in the thesis. The study benchmarks biological conservation, batch correction, and computational efficiency of various integration methods.

Below is a summary of the columns in the final integrated Atlas:

| Column                  | Type  | First 5 Values of the Columns                                                                                       |
|-------------------------|-------|-------------------------------------------------------------------------------------------------------------------|
| SampleID                | cells | OT_23-1_Pool_d20, OT_23-1_Pool_d20, OT_23-1_Pool_d20, OT_23-1_Pool_d20, OT_23-1_Pool_d20                           |
| ExperimentID            | cells | NI23.1, NI23.1, NI23.1, NI23.1, NI23.1                                                                             |
| donor                   | cells | HPSI0316i-aask_4, HPSI0316i-aask_4, HPSI0516i-oadp_4, HPSI0316i-ierp_4, HPSI0316i-ierp_4                           |
| Condition               | cells | Kabuki, Kabuki, Kabuki, Kabuki, Kabuki                                                                             |
| robustID                | cells | AAACCCAAGCAGCCTC_OT_23-1_Pool_d20, AAACCCAAGCGCAATG_OT_23-1_Pool_d20, AAACCCAGTCTTACTT_OT_23-1_Pool_d20, AAACCCATCAACGCTA_OT_23-1_Pool_d20, AAACCCATCATTCCTA_OT_23-1_Pool_d20 |
| original_labels         | cells | oRG, vRG, panRG-O, vRG, oRG                                                                                        |
| Day_fixed               | cells | 20.0, 20.0, 20.0, 20.0, 20.0                                                                                       |
| level_1                 | cells | neural_progenitor_cell, neural_crest, neural_progenitor_cell, neural_progenitor_cell, neural_crest                 |
| level_2                 | cells | mesencephalic_npc, nan, mesencephalic_npc, mesencephalic_npc, nan                                                  |
| level_3                 | cells | nan, nan, nan, nan, nan                                                                                            |
| level_4                 | cells | nan, nan, nan, nan, nan                                                                                            |
| level_5                 | cells | nan, nan, nan, nan, nan                                                                                            |
| final_level             | cells | mesencephalic_npc, neural_crest, mesencephalic_npc, mesencephalic_npc, neural_crest                               |
| Disease                 | cells | True, True, True, True, True                                                                                       |
| assay_sc                | cells | 10x 3' v3, 10x 3' v3, 10x 3' v3, 10x 3' v3, 10x 3' v3                                                              |
| batch                   | cells | 2d, 2d, 2d, 2d, 2d                                                                                                 |
| _scvi_batch             | cells | 0, 0, 0, 0, 0                                                                                                      |
| _scvi_labels            | cells | 10, 14, 10, 10, 14                                                                                                 |
| leiden                  | cells | 31, 31, 31, 16, 31                                                                                                 |
| conditions_combined     | cells | 2d, 2d, 2d, 2d, 2d                                                                                                 |
| highly_variable         | genes | True, True, True, True, True                                                                                       |
| means                   | genes | 2.6507799937333028, 1.0712801599771362, 2.3090139281158506, 0.04923861683770433, 0.3004041662210974               |
| dispersions             | genes | 3.4963160357675633, 4.243031764250064, 4.485409726957123, 2.543828781386361, 2.8509761056276655                    |
| dispersions_norm        | genes | -0.3227612, 1.8771468, 1.2808453, -0.30524805, 0.35034424                                                          |
| highly_variable_nbatches| genes | 2, 3, 3, 2, 2                                                                                                      |
| highly_variable_intersection | genes | False, True, True, False, False 
