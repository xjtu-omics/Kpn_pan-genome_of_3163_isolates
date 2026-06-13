from xgboost import XGBClassifier
from sklearn import metrics
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import metrics
import seaborn as sns

def model_train(train_data,test_data,save_path,med,type):

    # Features and target variable
    X = train_data.drop(columns=["label"])  # Features
    y = train_data["label"]  # Target variable
    X_test = test_data.drop(columns=["label"])  # Features
    y_test = test_data["label"]  # Target variable

    params = {
        'random_state': 42,
        'n_estimators': 400,
        'learning_rate': 0.05,
        'max_depth': 5,
        'subsample': 0.8,
        'colsample_bytree': 0.4,
        'n_jobs': -1,
    }

    model = XGBClassifier(**params)
    model.fit(X, y)

    # Extract feature importance
    feature_importances = pd.DataFrame({
        "Feature": X.columns,
        "Importance":  model.feature_importances_
    }).sort_values(by="Importance", ascending=False)
    #feature_importances = feature_importances[feature_importances['Importance'] != 0]
    feature_importances.to_csv(save_path+type+'_xgb_importance.csv', index=False, encoding='utf-8-sig')

    # Output model performance metrics
    y_pred=model.predict(X_test)
    with open(save_path+type+'_xgb_metrics.txt','w') as f:
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
    plt.savefig(save_path+type+'_xgb_conf_matrix.png')