# no_votes

import pandas as pd
import numpy as np


def foo(x):
    n_no_votes = len([e for e in x if e == 'NO_VOTE'])
    if n_no_votes > (len(x) / 2):
        return "* NO_VOTE *"
    else:
        return np.mean([e for e in x if e != 'NO_VOTE'])


data_arr = np.random.rand(100, 50)
data = pd.DataFrame(data_arr)

num_no_votes = 350

# Create Random Mask
rand_zero_one_mask = np.random.randint(2, size=data.shape)
# Fill df with 0 where mask is 0
data = data.where(rand_zero_one_mask == 0, "NO_VOTE")

data["NUM_NO_VOTES"] = data.apply(
    lambda x: len([e for e in x if e == 'NO_VOTE']), axis='columns'
)

data['VOTED_PREDICTION'] = data.apply(foo, axis=1)

print(data)
