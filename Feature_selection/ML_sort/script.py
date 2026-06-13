import pandas as pd
from sklearn.model_selection import train_test_split
import os
import RF,LR,SVM,XGboost
import numpy as np
def merge_label_to_fm(fm,label):

    fm = fm.copy()
    fm["label"] = label

    fm = fm[fm["label"].isin(["R", "S"])]
    fm["label"] = fm["label"].map({"R": 1, "S": 0})

    return fm

def tt_split(df,result_path,med,type):

    # Split features and labels
    X = df.drop('label', axis=1)  # Features
    y = df['label']  # Labels

    # Split at a 7:3 ratio while preserving class proportions
    X_train, X_test, y_train, y_test = train_test_split(
        X, y,
        test_size=0.3,
        random_state=42,  # Random seed for reproducibility
        stratify=y  # Stratified split by label
    )

    # Merge features and labels, then export to file
    train_set = pd.concat([X_train, y_train], axis=1)
    test_set = pd.concat([X_test, y_test], axis=1)

    if not os.path.exists(result_path+med+'/'):
        os.makedirs(result_path+med+'/')
    train_set.to_csv(result_path+med+'/'+type+'_train_set.csv')
    test_set.to_csv(result_path+med+'/'+type+'_test_set.csv')

    return train_set,test_set

if __name__ == '__main__':

    # Read R/S resistance phenotype data
    pheno_path="./Panaroo-DownStream-both/phenotypes.csv"
    phenos=pd.read_csv(pheno_path,index_col=0)
    # Read the names of 35 antibiotics
    #med_list=phenos.columns.tolist()
    med_list=['CFP','CMZ','CSL','CZA','DOX','MOX','NET','PIP','TMP']

    # Set the root directory for storing results
    result_path = './Panaroo-DownStream-both/ml_sort/'

    # Read two base feature matrices and put them in a list
    core_path = "./Panaroo-DownStream-both/final_core_feature_matrix.csv"
    core = pd.read_csv(core_path,index_col=0)
    dis_path = "./Panaroo-DownStream-both/final_dispensable_feature_matrix.csv"
    dis = pd.read_csv(dis_path,index_col=0)
    fm_dict = {'core':core,'dis':dis}

    for med in med_list:

        label=phenos[med]
        print(f"Processing : {med}")

        for type in fm_dict:

            # Merge labels into the feature matrix
            whole_data = merge_label_to_fm(fm_dict[type], label)
            # If split dataset result files already exist, skip the split and save steps.
            # Split X1 (train + test) and the validation set; save validation directly and use data for later ranking
            y = whole_data['label']
            mask_gn = whole_data.index.str.startswith("GN") # All IDs starting with GN
            # Anhui samples
            X_gn = whole_data[mask_gn]
            y_gn = y[mask_gn]
            # PATRIC samples
            X_non_gn = whole_data[~mask_gn]
            y_non_gn = y[~mask_gn]
            # Split PATRIC samples at 55%.
            try:
                if len(y_non_gn) == 0: #No PATRIC samples; validation is empty
                    X_non_gn_1 = pd.DataFrame(columns=X_non_gn.columns)
                    X_non_gn_2 = pd.DataFrame(columns=X_non_gn.columns)
                    y_non_gn_1 = pd.Series(dtype=y_non_gn.dtype)
                    y_non_gn_2 = pd.Series(dtype=y_non_gn.dtype)
                else:
                    unique, counts = np.unique(y_non_gn, return_counts=True)
                    stratify_y = y_non_gn if counts.min() >= 2 else None #If any resistance phenotype class has fewer than two samples, do not perform stratified sampling
                    X_non_gn_1, X_non_gn_2, y_non_gn_1, y_non_gn_2 = train_test_split(
                        X_non_gn, y_non_gn,
                        test_size=0.5,
                        random_state=42,
                        stratify=stratify_y
                    )

            except ValueError as e:
                # Other ValueError: re-raise or skip
                print(f"[split error] {e}")
                continue

            X1 = pd.concat([X_gn, X_non_gn_1]).astype(int)
            y1 = pd.concat([y_gn, y_non_gn_1]).astype(int)
            validation_fm=pd.concat([X_non_gn_2, y_non_gn_2], axis=1)

            # Split the training and test sets from X1
            train_fm,test_fm=tt_split(X1,result_path,med,type)
            y_train = train_fm['label']
            # After splitting, the training set may contain only one phenotype; in that case, discard this antibiotic.
            if y_train.nunique() < 2:
                print("Skip: only one class in y_train")
                continue

            #debug
            print(y_train.head())
            print(y_train.dtype)
            print(y_train.nunique())
            print(y_train.isna().sum())

            # Random forest model ranking
            RF.model_train(train_fm,test_fm,result_path+med+'/',med,type)
            # Logistic regression model ranking
            LR.model_train(train_fm,test_fm,result_path+med+'/',med,type)
            # LightGBM model ranking
            XGboost.model_train(train_fm,test_fm,result_path+med+'/',med,type)
            # SVM model ranking
            SVM.model_train(train_fm,test_fm,result_path+med+'/',med,type)

            # Save the validation set
            validation_fm.to_csv(result_path+med+'/'+type+'_validation_set.csv')