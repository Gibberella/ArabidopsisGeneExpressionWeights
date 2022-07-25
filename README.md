# ArabidopsisGeneExpressionWeights

**File Descriptions**

**model.xml**: the multi-tissue model of Arabidopsis thaliana's metabolic network used in the study

**WeightGeneration.ipynb**: a Jupyter Notebook that takes a model and omic dataset and generates weights for use in an FBA optimization

**FBA_and_FVA_Optimization.m**: a file that performs the FBA optimization that incorporates omic weights

**Example_Omic_Data.xlsx**: an example of the tissue-specific omic data format needed for the **WeightGeneration.ipynb** script

**Example_Missing_Data.xlsx**: an example of the missing omic data format needed for the **WeightGeneration.ipynb** script

**M_Protein_Weights.xlsx**: the weight file that was used, together with **model.xml**, to generate the protein-informed FBA solutions in the study
