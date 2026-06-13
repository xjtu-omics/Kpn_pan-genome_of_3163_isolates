from collections import Counter
import pandas as pd
import os

def process_ml_sort(med,core_dis,thresh):
    # For the specified feature-matrix type (core/dis) of a specified antibiotic, select features appearing in at least two of four train-set model rankings
    path=f'./Panaroo-DownStream-both/ml_sort/{med}/'

    if not os.path.exists(f"{path}{core_dis}_xgb_importance.csv"):
        return []

    # Read the data file
    all_elements=[]
    models=['xgb','lr','rf','svm']
    for model in models:
        df = pd.read_csv(f"{path}{core_dis}_{model}_importance.csv")
        # Get the first column
        first_column = df.iloc[:, 0]
        # Ensure k does not exceed the column length
        actual_t = min(thresh, len(first_column))
        # Get the first k elements and convert them to a list
        result_list = first_column.head(actual_t).tolist()
        # Append to the combined list
        all_elements.extend(result_list)

    # Count occurrences of each element
    element_counts = Counter(all_elements)
    print(element_counts)
    # Keep elements with occurrence count >= 2
    result = [elem for elem, count in element_counts.items() if count >= 2]
    print(result)

    return result

def build_fm(feature_list,fm):
    # Extract columns for each feature in feature_list plus the target column from fm, then assemble a smaller feature matrix
    new_fm = fm[feature_list + ['label']]
    new_fm = new_fm.astype(int)
    return new_fm