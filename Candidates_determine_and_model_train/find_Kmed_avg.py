###Final code file used to compute the results
#Compute the final K_med and corresponding feature set for each antibiotic
import pandas as pd
import json
import list2fm as lf
import ensemble_avg as en
import matplotlib.pyplot as plt
import os
import numpy as np
os.environ["PYTHONHASHSEED"] = "0"
os.environ["OMP_NUM_THREADS"] = "1"

# Read R/S resistance phenotype data
phenos=pd.read_csv("./Panaroo-DownStream-both/phenotypes.csv",index_col=0)

# Read antibiotic names
#med_list=phenos.columns.tolist()
med_list=['TZP']
K_list=[35]

# Read statistical-test filtering results
chi_kw_path="./Panaroo-DownStream-both/statistics-test/"
with open(chi_kw_path+'chi_and_kw_core.json', 'r', encoding='utf-8') as file:
    chi_and_kw_core = json.load(file)
with open(chi_kw_path+'chi_and_kw_dispensable.json', 'r', encoding='utf-8') as file:
    chi_and_kw_dis = json.load(file)

if __name__=='__main__':

    root_save_path="./Panaroo-DownStream-both/feature_set_determine/5f_avgauc_train/"
    root_read_path="./Panaroo-DownStream-both/"

    for i in range(len(med_list)):

        med = med_list[i]
        K_med = K_list[i]
        step_size = range(1, K_med + 1, 1)
        print(med)

        os.makedirs(f"{root_save_path}result/{med}/", exist_ok=True)

        # AUC result list
        auc_array=[]
        # Initialize feature_num before the loop
        last_total_and_fn=0
        # feature_num array
        feature_num=[]
        # Thresholds corresponding to feature_num
        fn_thresh=[]

        ### Read the core and dispensable train and test datasets for the current antibiotic
        whole_fm_path=f"{root_read_path}ml_sort/{med}/"
        train_core_fm=pd.read_csv(whole_fm_path+'core_train_set.csv',index_col=0)
        train_dis_fm=pd.read_csv(whole_fm_path+'dis_train_set.csv',index_col=0)
        test_core_fm=pd.read_csv(whole_fm_path+'core_test_set.csv',index_col=0)
        test_dis_fm=pd.read_csv(whole_fm_path+'dis_test_set.csv',index_col=0)

        for threshold in list(step_size):

            ### Combine statistical-test results with machine-learning rankings to filter features and obtain a feature-name list
            # ML filtering note: ranking results only retain weighted features. Check the relationship between the threshold and the total feature count.
            ml_core_list=lf.process_ml_sort(med,'core',threshold)
            ml_dis_list=lf.process_ml_sort(med,'dis',threshold)
            filted_core_list=list(set(ml_core_list) & set(chi_and_kw_core[med]))
            filted_dis_list=list(set(ml_dis_list) & set(chi_and_kw_dis[med]))
            filted_total_list=filted_core_list+filted_dis_list

            ### Use the feature-name list and original feature matrix to extract features and assemble a new feature matrix
            # Note: build_fm returns a feature matrix with the label column
            # Assemble the train matrix
            selected_train_core_fm=lf.build_fm(filted_core_list,train_core_fm)
            selected_train_dispensable_fm=lf.build_fm(filted_dis_list,train_dis_fm)
            selected_train_total_fm = pd.concat([selected_train_core_fm, selected_train_dispensable_fm], axis=1) # Merge
            selected_train_total_fm = selected_train_total_fm.loc[:, ~selected_train_total_fm.columns.duplicated()] # Deduplicate columns
            # Assemble the test matrix
            selected_test_core_fm=lf.build_fm(filted_core_list,test_core_fm)
            selected_test_dispensable_fm=lf.build_fm(filted_dis_list,test_dis_fm)
            selected_test_total_fm = pd.concat([selected_test_core_fm, selected_test_dispensable_fm], axis=1) # Merge
            selected_test_total_fm = selected_test_total_fm.loc[:, ~selected_test_total_fm.columns.duplicated()] # Deduplicate columns

            ### Five-fold cross-validation
            # Train only when both conditions are met: the feature matrix is not empty and the feature set has expanded compared with the previous run
            # After training, append new values to each result list and the corresponding fn_list
            # Update last_fn at the end
            if selected_train_total_fm.shape[1] > 1 and len(filted_total_list) > last_total_and_fn:
                total_auc=en.train_ensemble_with_cv(selected_train_total_fm,selected_test_total_fm,med)
                auc_array.append(total_auc)
                feature_num.append(len(filted_total_list))
                fn_thresh.append(threshold)
                last_total_and_fn=len(filted_total_list)

        if len(auc_array)==0:
            print(f"No valid feature set found for {med}. Skipping.")
            continue
        ### Find the maximum AUC value in auc_array and record it in a dictionary
        max_value = max(auc_array)  # Find the maximum value
        max_index = auc_array.index(max_value)  # Find the index of the first occurrence of the maximum value

        ### Plot
        # Create the canvas
        fig, ax = plt.subplots(figsize=(10, 6))

        # 1. Draw the line plot
        ax.plot(feature_num, auc_array, markersize=8)
        # Set integer ticks; adjust as needed for the specific case
        numpy_feature_num = np.array(feature_num)
        x_min, x_max = int(np.floor(numpy_feature_num.min())), int(np.ceil(numpy_feature_num.max()))
        plt.xticks(np.arange(x_min, x_max + 5, 5))

        # 2. Annotate the point and add dashed lines
        # Mark the point
        ax.scatter(feature_num[max_index], max_value, color='red', s=100, zorder=5)
        # Add dashed lines: a vertical x-axis line and a horizontal y-axis line
        ax.axvline(x=feature_num[max_index], color='gray', linestyle='--', alpha=0.5)
        ax.axhline(y=max_value, color='gray', linestyle='--', alpha=0.5)
        # Annotate values on the axes
        ax.text(feature_num[max_index], ax.get_ylim()[0]-0.004, f'x={feature_num[max_index]}', ha='center', va='top', color='blue')  # x-axis annotation
        ax.text(ax.get_xlim()[0]-1,max_value+0.005, f'y={max_value:.4f}', ha='right', va='center', color='green')  # y-axis annotation

        # 3. Polish the plot
        ax.set_xlabel('feature number')
        ax.set_ylabel('roc_auc value')
        ax.set_title('best cut-off point')
        ax.grid(True, linestyle=':', alpha=0.7)
        plt.tight_layout()
        plt.savefig(f'{root_save_path}result/{med}/bestkey.png')
        fig.clear()

        # Print the feature sets corresponding to all maxima
        # Output every maximum value
        # Write all results to one file
        ml_core_list=lf.process_ml_sort(med,'core',fn_thresh[max_index])
        ml_dis_list=lf.process_ml_sort(med,'dis',fn_thresh[max_index])
        feature_list=list(set(ml_core_list) & set(chi_and_kw_core[med])) + list(set(ml_dis_list) & set(chi_and_kw_dis[med]))
        with open(f'{root_save_path}result/{med}/feature_select_result.txt', 'w') as f:

            f.write("\n")
            f.write(f"roc_auc：{max_value}\n") # Corresponding AUC
            rank_path=f"{root_read_path}ml_sort/"
            xgb_core_rank=pd.read_csv(f'{rank_path}{med}/core_xgb_importance.csv')
            lr_core_rank=pd.read_csv(f'{rank_path}{med}/core_lr_importance.csv')
            rf_core_rank=pd.read_csv(f'{rank_path}{med}/core_rf_importance.csv')
            svm_core_rank=pd.read_csv(f'{rank_path}{med}/core_svm_importance.csv')
            xgb_dis_rank=pd.read_csv(f'{rank_path}{med}/dis_xgb_importance.csv')
            lr_dis_rank=pd.read_csv(f'{rank_path}{med}/dis_lr_importance.csv')
            rf_dis_rank=pd.read_csv(f'{rank_path}{med}/dis_rf_importance.csv')
            svm_dis_rank=pd.read_csv(f'{rank_path}{med}/dis_svm_importance.csv')

            core_chi_path=f"{root_read_path}statistics-test/core_chisquare.csv"
            core_chi=pd.read_csv(core_chi_path,index_col=0)
            core_kw_path=f"{root_read_path}statistics-test/core_kw.csv"
            core_kw=pd.read_csv(core_kw_path,index_col=0)
            dis_chi_path=f"{root_read_path}statistics-test/dispensable_chisquare.csv"
            dis_chi=pd.read_csv(dis_chi_path,index_col=0)
            dis_kw_path=f"{root_read_path}statistics-test/dispensable_kw.csv"
            dis_kw=pd.read_csv(dis_kw_path,index_col=0)

            for i in feature_list:
                print(i)
                if '->' in i:
                    #Core-gene features
                    chi_pvalue=core_chi.at[i,med]
                    kw_pvalue=core_kw.at[i,med]
                    xgb_weight=xgb_core_rank.loc[xgb_core_rank['Feature'] == i, 'Importance'].values[0]
                    if xgb_weight!=0:
                        xgb_rank=int(((xgb_core_rank.index[xgb_core_rank['Feature'] == i]).to_list())[0])+1
                    else:
                        xgb_rank=-1
                    lr_weight=lr_core_rank.loc[lr_core_rank['Feature'] == i, 'Importance'].values[0]
                    if lr_weight!=0:
                        lr_rank=int(((lr_core_rank.index[lr_core_rank['Feature'] == i]).to_list())[0])+1
                    else:
                        lr_rank=-1
                    rf_weight=rf_core_rank.loc[rf_core_rank['Feature'] == i, 'Importance'].values[0] # Already transformed with -log
                    if not pd.isna(rf_weight):
                        rf_rank=int(((rf_core_rank.index[rf_core_rank['Feature'] == i]).to_list())[0])+1
                    else:
                        rf_rank=-1
                    svm_weight=svm_core_rank.loc[svm_core_rank['Feature'] == i, 'Importance'].values[0]
                    if svm_weight!=0:
                        svm_rank=int(((svm_core_rank.index[svm_core_rank['Feature'] == i]).to_list())[0])+1
                    else:
                        svm_rank=-1
                else:
                    #Dispensable-gene features
                    chi_pvalue=dis_chi.at[i,med]
                    kw_pvalue=dis_kw.at[i,med]
                    xgb_weight=xgb_dis_rank.loc[xgb_dis_rank['Feature'] == i, 'Importance'].values[0]
                    if xgb_weight!=0:
                        xgb_rank=int(((xgb_dis_rank.index[xgb_dis_rank['Feature'] == i]).to_list())[0])+1
                    else:
                        xgb_rank=-1
                    lr_weight=lr_dis_rank.loc[lr_dis_rank['Feature'] == i, 'Importance'].values[0]
                    if lr_weight!=0:
                        lr_rank=int(((lr_dis_rank.index[lr_dis_rank['Feature'] == i]).to_list())[0])+1
                    else:
                        lr_rank=-1
                    rf_weight=rf_dis_rank.loc[rf_dis_rank['Feature'] == i, 'Importance'].values[0] # Already transformed with -log
                    if not pd.isna(rf_weight):
                        rf_rank=int(((rf_dis_rank.index[rf_dis_rank['Feature'] == i]).to_list())[0])+1
                    else:
                        rf_rank=-1
                    svm_weight=svm_dis_rank.loc[svm_dis_rank['Feature'] == i, 'Importance'].values[0]
                    if svm_weight!=0:
                        svm_rank=int(((svm_dis_rank.index[svm_dis_rank['Feature'] == i]).to_list())[0])+1
                    else:
                        svm_rank=-1

                f.write(f"{i},chi_pvalue={chi_pvalue},kw_pvalue={kw_pvalue},xgb_rank={xgb_rank},lr_rank={lr_rank},rf_rank={rf_rank},svm_rank={svm_rank}\n")

        ### For the final optimal feature set, call en.train_ensemble_with_cv again and save the model, related parameters, and ROC curve
        ## First assemble the final training feature matrix
        final_core_list=list(set(ml_core_list) & set(chi_and_kw_core[med]))
        final_dis_list=list(set(ml_dis_list) & set(chi_and_kw_dis[med]))
        # Assemble the train matrix
        final_train_core_fm=train_core_fm[final_core_list + ['label']]
        final_train_dis_fm=train_dis_fm[final_dis_list + ['label']]
        final_train_fm = pd.concat([final_train_core_fm, final_train_dis_fm], axis=1) # Merge
        final_train_fm = final_train_fm.loc[:, ~final_train_fm.columns.duplicated()] # Deduplicate columns
        # Assemble the test matrix
        final_test_core_fm=test_core_fm[final_core_list + ['label']]
        final_test_dis_fm=test_dis_fm[final_dis_list + ['label']]
        final_test_fm = pd.concat([final_test_core_fm, final_test_dis_fm], axis=1) # Merge
        final_test_fm = final_test_fm.loc[:, ~final_test_fm.columns.duplicated()] # Deduplicate columns
        # Call the function
        final_auc=en.train_ensemble_with_cv(final_train_fm,final_test_fm,med,model_save=True)
        print(f"Final AUC for {med} with K_med={K_med}: {final_auc}")
        # Save the train and test matrices
        final_train_core_fm.to_csv(f"{root_save_path}core_train_matrix/{med}.csv")
        final_train_dis_fm.to_csv(f"{root_save_path}dis_train_matrix/{med}.csv")
        final_test_core_fm.to_csv(f"{root_save_path}core_test_matrix/{med}.csv")
        final_test_dis_fm.to_csv(f"{root_save_path}dis_test_matrix/{med}.csv")
        final_train_fm.to_csv(f"{root_save_path}result_train_matrix/{med}.csv")
        final_test_fm.to_csv(f"{root_save_path}result_test_matrix/{med}.csv")