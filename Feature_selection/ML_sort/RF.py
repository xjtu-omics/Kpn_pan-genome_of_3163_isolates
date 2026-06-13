import pandas as pd
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt
from sklearn import metrics
import seaborn as sns
import math

def model_train(train_data,test_data,save_path,med,type):

    # Features and target variable
    X = train_data.drop(columns=["label"])  # Features
    y = train_data["label"]  # Target variable
    X_test = test_data.drop(columns=["label"])
    y_test = test_data["label"]

    rf_model = RandomForestClassifier(random_state=42, max_features='sqrt', n_estimators=500, n_jobs=-1)
    rf_model.fit(X, y)  # Train using the full dataset

    # 3. Extract feature importance and sort by ascending -log10 p-value
    feature_importances = pd.DataFrame({
        "Feature": X.columns,
        "Importance":  rf_model.feature_importances_
    }).sort_values(by="Importance", ascending=False)
    for index,row in feature_importances.iterrows():
        #x=float(math.log10(row["Importance"]))
        if row["Importance"]>0:
            feature_importances.at[index,"Importance"]=-math.log10(row["Importance"])
        else:
            feature_importances.at[index,"Importance"]=float('nan')

    feature_importances.to_csv(save_path+type+'_rf_importance.csv', index=False, encoding='utf-8-sig')

    # Output model performance metrics
    y_pred=rf_model.predict(X_test)
    with open(save_path+type+'_rf_metrics.txt','w') as f:
        f.write('accuracy: '+str(metrics.accuracy_score(y_test,y_pred))+'\n')
        f.write('precision: '+str(metrics.precision_score(y_test,y_pred))+'\n')
        f.write('recall: '+str(metrics.recall_score(y_test,y_pred))+'\n')
        f.write('f1: '+str(metrics.f1_score(y_test,y_pred))+'\n')

    # Output the confusion matrix
    conf_matrix = metrics.confusion_matrix(y_test, y_pred, labels=[0,1])
    # Visualize
    plt.figure(figsize=(8, 6))
    sns.heatmap(conf_matrix, annot=True, fmt='d', cmap='Blues',
                xticklabels=['Predicted Sensitive', 'Predicted Resistant'],
                yticklabels=['Actual Sensitive', 'Actual Resistant'])
    plt.xlabel('Predicted Label')
    plt.ylabel('True Label')
    plt.title(med+' Confusion Matrix')
    plt.savefig(save_path+type+'_rf_conf_matrix.png')