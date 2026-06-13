import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['font.family'] = 'Arial'

def plot_metrics(data):

    # Prepare data
    keys = list(data.keys())
    means = [np.mean(data[k]) for k in keys]
    print(means)
    stds = [np.std(data[k], ddof=1) for k in keys]  # Sample standard deviation

    colors = [
        # CZO, FOX, CXM, CAZ
        '#377EB8', '#377EB8', '#377EB8', '#377EB8',
        # ETP, IPM, MEM
        '#B23648', '#B23648', '#B23648',
        # SAM, TZP
        '#DC7369', '#DC7369',
        # GEN, AMK, TOB
        '#D8EBCD', '#D8EBCD', '#D8EBCD',
        # TCY
        '#F8EFB5',
        # CIP, LVX
        '#DAD4B9', '#DAD4B9',
        # NIT
        '#C8CDCF',
        # SXT
        '#E1F3FA'
    ]

    # x-axis indices
    x = np.arange(len(keys))

    # Create the plot
    plt.figure(figsize=(12, 6))

    # Draw the bar plot
    plt.bar(x, means, yerr=stds, capsize=5, color=colors)

    # Set axes
    plt.xticks(x, keys)
    plt.xlabel('Antibiotics', fontsize=14)
    plt.ylabel('roc-auc', fontsize=14)
    plt.title('Test', fontsize=16)

    # Modify the outer spine width
    plt.gca().spines['top'].set_linewidth(1.3)
    plt.gca().spines['right'].set_linewidth(1.3)
    plt.gca().spines['bottom'].set_linewidth(1.3)
    plt.gca().spines['left'].set_linewidth(1.3)

    # Add a gray dashed line at y=0.8
    plt.axhline(y=0.8, color='gray', linestyle='--')

    # Force the y-axis upper bound to 1
    plt.ylim(bottom=0, top=1)

    # Adjust layout
    plt.tight_layout()

    # Show the plot
    # 2) Change the file name here to include the metric type and dataset type
    plt.savefig("./Figures-re/fig4/4b-test.pdf")


if __name__ == "__main__":

    test_auc={'CZO':0.94, 'FOX':0.87, 'CXM':0.96, 'CAZ':0.95, 'ETP':0.86, 'IPM':0.94, 'MEM':0.91,
              'SAM':0.88, 'TZP':0.90, 'GEN':0.93, 'AMK':0.88, 'TOB':0.93, 'TCY':0.92, 'CIP':0.94,
              'LVX':0.96, 'NIT':0.87, 'SXT':0.93}
    '''
    validation_auc={'CZO':0.96, 'FOX':0.92, 'CXM':0.98, 'CAZ':0.92, 'ETP':0.91, 'IPM':0.97, 'MEM':0.95,
                    'SAM':0.93, 'TZP':0.83, 'GEN':0.97, 'AMK':0.91, 'TOB':0.98, 'TCY':0.91, 'CIP':0.99,
                    'LVX':0.97, 'NIT':0.80, 'SXT':0.94}
    '''

    plot_metrics(test_auc)