import numpy as np


def cal_prior(adata_train_dummy, cdata_train_dummy, a_label, c_label, dat_num, m, p):

    prior_probability = np.zeros((len(c_label), len(a_label)))


    for j in range(0, len(c_label)):
        denominator = 0
        for i in range(0, dat_num):
            if cdata_train_dummy[i, j] == 1:
                denominator += 1
                prior_probability[j, :] += adata_train_dummy[i, :]
            else:
                continue;
        prior_probability[j, :] += m*p
        prior_probability[j, :] /= (denominator + m)

    return  prior_probability


def cal_prob(adata_train_dummy, dat_num):
    probability = np.sum(adata_train_dummy, axis=0)
    probability = probability/dat_num
    return probability