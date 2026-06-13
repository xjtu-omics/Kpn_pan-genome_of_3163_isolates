import pandas as pd
from scipy.stats import kruskal
from toolkit import merge_label_to_fm

# Read R/S resistance phenotype data
pheno_path="./Panaroo-DownStream-both/phenotypes.csv"
phenos=pd.read_csv(pheno_path,index_col=0)

fm_path="./Panaroo-DownStream-both/final_core_feature_matrix.csv"
fm=pd.read_csv(fm_path,index_col=0)

# Remove zero-variance features where all values are identical
fm = fm.loc[:, fm.nunique() > 1]

med_list=phenos.columns.tolist()
variants=fm.columns.tolist()
result=pd.DataFrame(columns=med_list,index=variants)

# Get labels
#med='AMP'
for med in med_list:

    RS_data=phenos[med]
    # Merge labels with the large variant matrix
    df=merge_label_to_fm(fm,RS_data)

    # Assume df is a DataFrame containing the feature matrix and y is the target variable
    X = df.drop(columns='label')  # Delete the target-variable column to obtain the feature matrix
    y = df['label']  # Target variable

    # Remove features whose column values are all identical
    X = X.loc[:, X.nunique() > 1]  # Keep only features with more than one distinct value

    # Store the p-value for each feature
    p_values = []

    for feature in X.columns:

        # Group by phenotype
        grouped = [X[feature][y == c] for c in y.unique()]

        # Keep only groups with at least two samples and within-group variation
        valid_groups = [g for g in grouped if len(g) > 1 and g.nunique() > 1]

        if len(valid_groups) >= 2:
            # Run the Kruskal-Wallis test
            _, p_value = kruskal(*valid_groups)
        else:
            p_value = float('nan')

        p_values.append(p_value)
        result.loc[feature, med] = p_value

    # Combine p-values and feature names into a DataFrame
    p_values_df = pd.DataFrame({
        'Feature': X.columns,
        'P-Value': p_values
    })

result.to_csv("./Panaroo-DownStream-both/statistics-test/core_kw.csv")