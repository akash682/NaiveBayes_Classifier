import numpy as np
def original_max(sequence):
    maximum = sequence[0]
    for item in sequence:
        if item[1] > maximum[1]:
            maximum = item
    list_max = []
    for item in sequence:
        if item[1] == maximum[1]:
            list_max.append(item)
    return list_max

def naivebayes(adata_test_dummy, prior_probability, a_probability, c_probability, a_label, c_label):
    oracle = []
    for item in c_label:
        oracle.append([item, 0])

    for i in range(0, len(c_label)):
        posterior_probability = 1
        likelihood = 1
        for j in range(0, len(a_label)):
            if adata_test_dummy[0, j] == 1:
                posterior_probability *= prior_probability[i, j]
                likelihood *= a_probability[j]
            elif adata_test_dummy[0,j] == 0:
                continue

        oracle[i][1] = posterior_probability
    max_data = original_max(oracle)
    return max_data



