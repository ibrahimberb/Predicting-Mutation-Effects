# no_votes

import pandas as pd
import numpy as np


def foo(x):
    n_no_votes = len([e for e in x if e == 'NO_VOTE'])
    if n_no_votes >= (len(x) / 2):
        return "NO_VOTE"
    else:
        return round(np.mean([e for e in x if e != 'NO_VOTE']), 5)


def bar(x):
    if x == 'NO_VOTE':
        return 'NO_VOTE'

    return int(float(x) >= 0.50)


data_arr = np.random.rand(100, 4)
data = pd.DataFrame(data_arr)

num_no_votes = 350

# Create Random Mask
rand_zero_one_mask = np.random.randint(2, size=data.shape)
# Fill df with 0 where mask is 0
data = data.where(rand_zero_one_mask == 0, "NO_VOTE")

data['PROB_1s_AVG'] = data.apply(foo, axis='columns')

# data["NUM_NO_VOTES_(ensure)"] = data.apply(
#     lambda x: len([e for e in x if e == 'NO_VOTE']), axis='columns'
# )

data['V1'] = data['PROB_1s_AVG'].apply(
    lambda x: 'NO_VOTE' if x == 'NO_VOTE' else int(float(x) >= 0.50)
)

data['V2'] = data['PROB_1s_AVG'].apply(
    bar
)

print(data)
