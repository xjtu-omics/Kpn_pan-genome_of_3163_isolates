import pandas as pd
from sklearn.linear_model import LogisticRegression
import matplotlib.pyplot as plt
from sklearn import metrics
import seaborn as sns

def model_train(train_data,test_data,save_path,med,type):

    # Features and target variable
    X = train_data.drop(columns=["label"])  # Features
    y = train_data["label"]  # Target variable
    X_test = test_data.drop(columns=["label"])
    y_test = test_data["label"]

    # Create the logistic regression model
    model = LogisticRegression(random_state=42, penalty='l1', solver='liblinear', C=0.1, max_iter=5000, n_jobs=-1)

    # Train the model
    model.fit(X, y)

    # Predict on the test set
    y_pred = model.predict(X_test)
    # Write prediction results to file
    with open(save_path+type+'_lr_metrics.txt','w') as f:
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
    plt.savefig(save_path+type+'_lr_conf_matrix.png')

    # Output weights for each feature
    feature_weights = model.coef_
    feature_names = X.columns.tolist()  # Replace with actual feature names
    weights = feature_weights[0]  # Get the weight array
    # Save feature names and weights to a DataFrame
    df = pd.DataFrame({
        'Feature': feature_names,
        'Importance': weights
    })
    # Sort
    df = df.sort_values(by='Importance', key=lambda col: col.abs(), ascending=False)
    #df = df[~df['Importance'].isin([0, 0.0])]
    df.to_csv(save_path+type+'_lr_importance.csv', index=False, encoding='utf-8-sig')