import numpy as np
import preprocess_data as p
import readfile as r
import cal_probability as cd
import naivebayes as nb

data_train = r.readfile("Text_data_for_Project1_train_data.txt")
data_test = r.readfile("Text_Data_for_Project1_test_data.txt")

alabel_train, clabel_train, adata_train, cdata_train = p.p_traindata(data_train)
alabel_test, adata_test = p.p_testdata(data_test)

atraintestmix = np.vstack((adata_train, adata_test))

att_num = np.shape(adata_train)[1]
dat_num = np.shape(adata_train)[0]

adata_traintest_dummy, a_label = p.indexing(atraintestmix, att_num, dat_num, key="test")

adata_train_dummy = adata_traintest_dummy[:-1, :]
adata_test_dummy = np.array([adata_traintest_dummy[-1, :]])

cdata_train_dummy, c_label = p.indexing(cdata_train, 1, dat_num, key="test")

m, p = 0, 0
prior_probability = cd.cal_prior(adata_train_dummy, cdata_train_dummy, a_label, c_label, dat_num, m, p)

a_probability = cd.cal_prob(adata_train_dummy, dat_num)
c_probability = cd.cal_prob(cdata_train_dummy, dat_num)

result = nb.naivebayes(adata_test_dummy, prior_probability, a_probability, c_probability, a_label, c_label)
print(result)