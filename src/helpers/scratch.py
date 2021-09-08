# import pandas as pd
# import numpy as np
#
#
# def convert_prob_0_to_class(prob_0):
#     """
#     Converts the value of probability of 0 into corresponding class.
#     If prob_0 greater than or equals to 0.50, it implies that the class is 0.
#     Otherwise, class is 1.
#     """
#
#     if prob_0 >= 0.50:
#         return 0
#     else:
#         return 1
#
#
# n_rows, n_cols = 100, 10
#
# data_arr = np.random.rand(n_rows, n_cols)
# data = pd.DataFrame(data_arr)
#
# data.columns = [f'col_{i}' for i in range(n_cols)]
# data['V_PREDICTION'] = data['col_9'].apply(lambda x: convert_prob_0_to_class(x))
# data['V_PREDICTION_2'] = data['col_9'].apply(convert_prob_0_to_class)
# print(data['V_PREDICTION'].equals(data['V_PREDICTION_2']))
#
# print(data)

import pandas as pd
from src.helpers.helpers_predator.displayers import visualize_label_counts
import matplotlib.pyplot as plt

data = pd.DataFrame({
    'VAL': [1, 2, 3, 4, 5],
    'Mut_eff_label': [1, 0, 0, 0, 0]
})

print(data)
val_counts = data['Mut_eff_label'].value_counts().sort_index()
print(val_counts)
val_counts = val_counts.rename({0: 'Disrupting', 1: 'Increasing + No Effect'})
print(val_counts)
print(data['Mut_eff_label'].value_counts().sort_index())
print(data['Mut_eff_label'].value_counts().index)
visualize_label_counts(data, 'Mut_eff_label')
plt.show()

# from common import export_data
# import pandas as pd
# import numpy as np
#
# data = pd.DataFrame(np.random.randn(100, 5))
# print(data)
#
# export_data(data, 'filename')
#
