import pandas as pd
import numpy as np


def convert_prob_0_to_class(prob_0):
    """
    Converts the value of probability of 0 into corresponding class.
    If prob_0 greater than or equals to 0.50, it implies that the class is 0.
    Otherwise, class is 1.
    """

    if prob_0 >= 0.50:
        return 0
    else:
        return 1


n_rows, n_cols = 100, 10

data_arr = np.random.rand(n_rows, n_cols)
data = pd.DataFrame(data_arr)

data.columns = [f'col_{i}' for i in range(n_cols)]
data['V_PREDICTION'] = data['col_9'].apply(lambda x: convert_prob_0_to_class(x))
data['V_PREDICTION_2'] = data['col_9'].apply(convert_prob_0_to_class)
print(data['V_PREDICTION'].equals(data['V_PREDICTION_2']))

print(data)
