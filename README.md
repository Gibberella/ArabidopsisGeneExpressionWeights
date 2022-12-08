# ArabidopsisGeneExpressionWeights

**File Descriptions**

**model.xml**: the multi-tissue model of Arabidopsis thaliana's metabolic network used in the study

**WeightGeneration.ipynb**: a Jupyter Notebook that takes a model and omic dataset and generates weights for use in an FBA optimization

**FBA_and_FVA_Optimization.m**: a file that performs the FBA optimization that incorporates omic weights

**Example_Omic_Data.xlsx**: an example of the tissue-specific omic data format needed for the **WeightGeneration.ipynb** script

**Example_Missing_Data.xlsx**: an example of the missing omic data format needed for the **WeightGeneration.ipynb** script

**M_Protein_Weights.xlsx**: the weight file that was used, together with **model.xml**, to generate the protein-informed FBA solutions in the study

**Supplemental_Dataset_1.xlsx**: Mapping between FBA and 13C-MFA/kinetic flux profiling studies' reactions.

**Supplemental_Dataset_2.xls**: The multi-tissue model of Arabidopsis thaliana's metabolic network used in the study in XLS format.

**Supplemental_Dataset_3.xlsx**: All Flux Variability Analysis results from the present study.

**Supplemental_Dataset_4.xlsx**: All flux vectors generated from the present study.

**Supplemental_Dataset_5.xlsx**: All weighted average error calculations done for the present study.

**Supplemental_Dataset_6.xlsx**: All weight vectors generated and used in the present study.
