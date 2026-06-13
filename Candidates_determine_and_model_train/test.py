import os
import numpy as np
import pandas as pd
from sklearn.model_selection import cross_validate, StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, VotingClassifier
from sklearn.svm import SVC
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import confusion_matrix
from xgboost import XGBClassifier
from sklearn.metrics import make_scorer, f1_score, matthews_corrcoef, accuracy_score, precision_score, recall_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score
import joblib

def model_test(med, labeled_data):

    # Import the feature matrix used during training (original training data or column names)
    train_df = pd.read_csv(f"./Panaroo-DownStream-patric/feature_set_determine/5f_avgauc_train/result_train_matrix/{med}.csv", index_col=0)  # Feature matrix used during training
    train_df = train_df.sort_index(axis=1)
    model_features = list(train_df.columns)  # Keep the column order used during training

    # Import the trained model
    model = joblib.load(f"./Panaroo-DownStream-patric/feature_set_determine/5f_avgauc_train/result/{med}/model_train_final.pkl")

    # Predict
    X = labeled_data[model_features]
    X = X.drop(columns=["label"])
    y_true = labeled_data["label"]  # True label column
    y_pred = model.predict(X)

    # Compute the confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    print("Confusion matrix:")
    print(cm)

    # Some metrics require probabilities
    if hasattr(model, "predict_proba"):
        y_proba = model.predict_proba(X)[:, 1]
    else:
        # For models that do not support predict_proba, such as some SVMs, use decision_function
        y_proba = model.decision_function(X)

    # Compute metrics and store them in a dictionary
    metrics_dict = {}
    if y_true is not None and len(set(y_true)) > 1:
        metrics_dict = {
            "F1": f1_score(y_true, y_pred, average="weighted"),
            "MCC": matthews_corrcoef(y_true, y_pred),
            "AUC": roc_auc_score(y_true, y_proba),
            "Accuracy": accuracy_score(y_true, y_pred),
            "Precision": precision_score(y_true, y_pred, average="weighted"),
            "Recall": recall_score(y_true, y_pred, average="weighted"),
        }
        print(metrics_dict)
    elif len(set(y_true)) == 1:
        metrics_dict = {
            "F1": f1_score(y_true, y_pred, average="weighted"),
            "MCC": matthews_corrcoef(y_true, y_pred),
            "AUC": None,
            "Accuracy": accuracy_score(y_true, y_pred),
            "Precision": precision_score(y_true, y_pred, average="weighted"),
            "Recall": recall_score(y_true, y_pred, average="weighted"),
        }
    else:
        print("No label column; metrics cannot be computed.")

    return metrics_dict

if __name__ == "__main__":

    # Get the antibiotic list
    feature_folder_path="./Panaroo-DownStream-patric/feature_set_determine/5f_avgauc_train/result/"
    folder_list = [f for f in os.listdir(feature_folder_path)]
    #folder_list = ['ETP'] # Run NIT only

    # Use the feature matrix corresponding to the test set for training
    for med in folder_list:
        print(med)
        df=pd.read_csv(f"./Panaroo-DownStream-patric/feature_set_determine/5f_avgauc_train/result_test_matrix/{med}.csv", index_col=0)
        result=model_test(med, df)
        with open(f"./Panaroo-DownStream-patric/feature_set_determine/5f_avgauc_train/test_result/{med}.txt", "w") as f:
            for key, value in result.items():
                f.write(f'{key}: {value}\n')