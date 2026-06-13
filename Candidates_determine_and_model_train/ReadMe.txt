(1) Train set
    ↓
    5-fold CV
    ↓
    Select the model / hyperparameters

(2) Refit the best model on the full train set

(3) Independent test set
    ↓
    Single prediction
    ↓
    Single AUC + ROC

ensemble_avg.py:
Ensemble model training function
Saved files: 1) metrics dataframe, 2) predicted probabilities and labels for each sample, 3) ROC-AUC curve, 4) trained model (.pkl format)
All files saved directly by soft_ensemble must be called and saved explicitly at the end; otherwise, the final output will correspond to the last K value in the loop.
(Here, five-fold cross-validation is performed on 30% of the samples)

list2fm.py:
Function for assembling a feature matrix from a feature-name list

find_Kmed_avg.py
This code file should be run twice.
First run: set the interval for all 35 antibiotics to 100, inspect the K-AUC trend, and determine the interval where the maximum should be selected, assuming AUC reaches 0.8.
Second run: find the maximum within the selected interval.
Saved files: 1) feature list and maximum AUC value, 2) K-AUC curve

selected_fm.py:
Based on the selected features, output the feature matrix corresponding to each antibiotic significant feature set

model_save.py:
Use a specific feature matrix as input and output the trained model and corresponding metrics

test.py：
Evaluate the model trained on 70% of the samples using the 30% test set