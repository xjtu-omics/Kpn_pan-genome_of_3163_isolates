import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

# Data
data = np.array([
    [1768, 7],
    [119, 7]
])

# Create the figure with width greater than height for a flatter rectangle
fig, ax = plt.subplots(figsize=(8, 3))   # <<< Adjust the wide aspect ratio here

ax.axis('off')

# Create the table
table = ax.table(
    cellText=data,
    cellLoc='center',
    loc='center'
)

# Adjust cell size
table.scale(2.1, 2.5)   # <<< The first value scales width and the second scales height

# Polish the style
for (row, col), cell in table.get_celld().items():
    cell.set_edgecolor('black')
    cell.set_linewidth(2)
    cell.set_fontsize(20)

plt.tight_layout()

# To save as PDF:
save_path = "./Figures-re\\fig5"
plt.savefig(save_path + "\\5d-ompA-NIT-LorData.pdf", format="pdf", bbox_inches="tight")

plt.show()
