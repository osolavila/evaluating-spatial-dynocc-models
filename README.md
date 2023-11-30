# evaluating-spatial-dynocc-models
Analyses supporting the publication: 'Evaluating the influence of neighborhood connectivity and habitat effects in dynamic occupancy species distribution models' Oriol Solà, Núria Aquilué, Sara Fraixedas, Lluís Brotons.

Exploring new approaches and methodologies to characterize species distribution dynamics, instead of solely relying on static spatial patterns, should be a priority for species distribution modeling research. Dynamic occupancy models (here, 'dynocc models') are a promising tool to capture temporal patterns of distribution change but their spatial accuracy has been shown to be limited. In this study, we evaluated the effectiveness of incorporating neighborhood connectivity effects into the colonization and extinction functions of dynocc models. For that, we compared dynocc models accounting either for neighborhood connectivity only, for site-level habitat covariates only, or combining both neighborhood and habitat explanations in the models for species extinction and colonization. All models were evaluated for a total of 46 bird species typical of forests and shrublands using data at 1km2 scale from two Catalan breeding bird atlases (CBBA2: 1999-2002) and (CBBA3: 2015-2018). We found that for most bird species (80%), dynocc models with habitat covariates alone predicted better than models only accounting for connectivity. Combining connectivity and habitat covariates improved model performance even more, especially for species with similar performance gains for habitat and connectivity models alone. This study shows the benefits of considering more spatially explicit formulations in dynocc models, specifically incorporating neighborhood connectivity into the extinction and colonization processes. Our work also highlights the importance of evaluating different model formulations and assessing which aspects of the model are more important depending on the study species.


Script files to replicate the data preparation, model fitting, analysis, and plot generation are provided here. Species occurrence and environmental data used for fitting models are available upon reasonable request (contact: o.sola@creaf.uab.cat).



Metadata:
1 Models_CBBA_species
-	Models_CBBA_species/Functions: help functions to support model fitting and evaluation.
-	Models_CBBA_species/Data_preparation: loads CBBA2 and CBBA3 species occurrence data and supporting data for the analysis.
-	Models_CBBA_species/model_fitting: Scripts to fit 2 season dynocc models to CBBA2 and 3. To reproduce the analysis, scripts should be run in numerical order. “Functions” and “Data_preparation” scripts are automatically loaded. Models fits are outputted as “.Rdata” in the specified “output.path” variable.
-	Models_CBBA_species/plots_and_evaluation: scripts to reproduce figures 1, 2, S1, S2, S3 and tables 1, S2, S3 from the manuscript. Loads model fits “.Rdata” files and performs the model evaluation as described in the manuscript.

2 Models_simulation_study
-	01_simulate_occu.R: defines the study area and sampling protocol and simulates occurrence for the different colonization and extinction scenarios as describe in the manuscript.
-	02_fit_models.R: Fits models to the simulate occupancy dynamics data. Models fits are outputted as “.Rdata” in the specified “output.path” variable.
-	03_evaluate_models.R: script to reproduce figure 3 from the manuscript. Loads model fits “.Rdata” files and performs the model evaluation as described in the manuscript.

3 Colext_IFM_example
Code and worked example for the fitting function of the spatial IFM colext. It includes example species occurrence data.
