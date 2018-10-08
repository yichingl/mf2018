import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import numpy as np
import pandas as pd

from sklearn import preprocessing
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

# return columns to drop
def get_to_drop(data_dataframe, decorr_to_drop):

	# get correlation information
	corr = data_dataframe.corr(method='pearson').abs(); # matrix form
	# from: https://chrisalbon.com/machine_learning/feature_selection/drop_highly_correlated_features/
	upper = corr.where(np.triu(np.ones(corr.shape), k=1).astype(np.bool));
	to_drop = [column for column in upper.columns if any(upper[column] > 0.95)];

	with open(decorr_to_drop,'w') as f:
		for drop_name in to_drop:
			f.write(drop_name+'\n');
	
	return to_drop;




##############################
# GROUP AND CTL NAME INPUT
##############################
# Group name
group_name = 'allcat_bal';
# Define train and test names for each class
boat_train_ctl_name = 'allcat_bal.boat';
heli_train_ctl_name = 'allcat_bal.heli';
boat_test_ctl_name = 'allcat_bal.boat';
heli_test_ctl_name = 'allcat_bal.heli';
# For GMM-UBM
# for training background models
traingmm_ctl_names = ['allcat_bal.boat','allcat_bal.heli']; 
# create list of train/test files to use to extract supervectors
gensupervect_train_ctl_names = [boat_train_ctl_name, heli_train_ctl_name];
gensupervect_test_ctl_names = [boat_test_ctl_name, heli_test_ctl_name];

##############################
# TOGGLE TASKS
##############################
FIND_PCA_COMP = 0;
GEN_VIS_PCA = 1;

##############################
# GMM PARAMETERS
##############################
n_clusters = 16;

##############################
# SPECGRAM TYPE
##############################
feat_types = ['mod','spec','mel'];
feat_type = 'spec';

if feat_type == 'spec':
	n_feats = 133;
	n_timeframes = 669;

	# [ 71.55  88.96  94.46  96.5   97.57  98.47  98.93  99.29  99.52  99.68
    #   99.82  99.92  99.96  99.98  99.99 100.   100.01 100.01 100.01 100.01
 	#   100.01 100.01]
	n_PCA_comps = 15;

elif feat_type == 'mel':
	n_feats = 128;
	n_timeframes = 173;

elif feat_type == 'mod':
	n_feats = 128;
	n_timeframes = 320;
	

##############################
# FULL_NAME
##############################
full_name = group_name+'.'+feat_type+'.'+str(n_clusters)+'n_cl'+'.'+str(n_feats)+'n_feats';

# Writing supervect path
supervect_out = '/home/yichingl/mf/supervect_out/';
# Dataframe output path
train_dataframe_out_filepath = supervect_out+full_name+'.train_dataframe.csv';
test_dataframe_out_filepath = supervect_out+full_name+'.test_dataframe.csv';
# Uncorrelated Dataframe to drop output path
decorr_to_drop = supervect_out+full_name+'.decorr_to_drop.csv';
# Uncorrelated Dataframe output path
train_dataframe_decorr_out_filepath = supervect_out+full_name+'.train_dataframe.decorr.csv';
test_dataframe_decorr_out_filepath = supervect_out+full_name+'.test_dataframe.decorr.csv';
# Visualize PCA cumulative variance explained
train_cum_var_PCA_impath = '../vis/'+full_name+'.cumvarPCA.png';
# Dataframe dim reduced output path
train_dataframe_decorr_pca_filepath = supervect_out+full_name+'.train_dataframe.pca.csv';
test_dataframe_decorr_pca_filepath = supervect_out+full_name+'.test_dataframe.pca.csv';
# Save representation using TSNE to csv for plotting
train_tsne_3D_data = '../vis/'+full_name+'.TSNE.3D.csv';
train_tsne_2D_data = '../vis/'+full_name+'.TSNE.2D.csv';


##############################
# CODE
##############################
# Read Dataframe
train_dataframe = pd.read_csv(train_dataframe_out_filepath, index_col=0);
# Extract labels and subgroups
labels_train = train_dataframe['labels_train'];
subgroups = train_dataframe['subgroups'];
train_dataframe = train_dataframe.drop(['labels_train','subgroups'],axis=1);
# Get drop features
to_drop = get_to_drop(train_dataframe,decorr_to_drop);
# Drop features and save decorrelated data
train_decorr_dataframe =  train_dataframe.drop(to_drop,axis=1)
train_decorr_dataframe = train_decorr_dataframe.assign(labels_train=labels_train, subgroups=subgroups);
train_decorr_dataframe.to_csv(train_dataframe_decorr_out_filepath);

# Read test set and drop the same features
test_dataframe = pd.read_csv(test_dataframe_out_filepath, index_col=0);
# Extract labels and subgroups
labels_test = test_dataframe['labels_test'];
test_dataframe = test_dataframe.drop(['labels_test'],axis=1);
# Drop features and save decorrelated data
test_decorr_dataframe =  train_dataframe.drop(to_drop,axis=1)
test_decorr_dataframe = test_decorr_dataframe.assign(labels_test=labels_test);
test_decorr_dataframe.to_csv(test_dataframe_decorr_out_filepath);

# Fit to train data
train_decorr_dataframe = train_decorr_dataframe.drop(['labels_train','subgroups'],axis=1);
train_data = train_decorr_dataframe.values;
test_decorr_dataframe = test_decorr_dataframe.drop(['labels_test'],axis=1);
test_data = test_decorr_dataframe.values;

# Scale data
scaler = preprocessing.StandardScaler();
scaler.fit_transform(train_data);

# Perform PCA on train data
if FIND_PCA_COMP:
	pca = PCA();
	pca.fit_transform(train_data);

	# print performance
	var = pca.explained_variance_ratio_;
	var_cum = np.cumsum(np.round(var,decimals=4)*100);
	print(var_cum);

	fig = plt.figure()
	plt.plot(var_cum);
	plt.xlabel('Component')
	plt.ylabel('Cumulative Explained Variance')
	fig.savefig(train_cum_var_PCA_impath);
if GEN_VIS_PCA:
	# perform PCA on train and test data
	pca = PCA(n_components=n_PCA_comps);
	pca.fit_transform(train_data);
	pca.transform(test_data);

	# print performance
	var = pca.explained_variance_ratio_;
	var_cum = np.cumsum(np.round(var,decimals=4)*100);
	print(var_cum);

	# TRAIN DATA
	# save dim red data
	# generate dataframe row and col labels
	sample_inds = np.arange(np.shape(train_data)[0]);
	feat_names = ['comp_'+str(n) for n in np.arange(np.shape(train_data)[1])];
	# create dataframe object
	train_dataframe = pd.DataFrame(data=train_data,index=sample_inds, columns=feat_names);
	train_dataframe = train_dataframe.assign(labels_train=labels_train, subgroups=subgroups);
	train_dataframe.to_csv(train_dataframe_decorr_pca_filepath);

	# TEST DATA
	# save dim red data
	# generate dataframe row and col labels
	sample_inds = np.arange(np.shape(test_data)[0]);
	feat_names = ['comp_'+str(n) for n in np.arange(np.shape(test_data)[1])];
	# create dataframe object
	test_dataframe = pd.DataFrame(data=test_data,index=sample_inds, columns=feat_names);
	test_dataframe = test_dataframe.assign(labels_train=labels_test);
	test_dataframe.to_csv(test_dataframe_decorr_pca_filepath);


	feat_embedded = TSNE(n_components=3).fit_transform(train_data)
	feat_embedded = np.asarray(feat_embedded)
	np.savetxt(train_tsne_3D_data, feat_embedded, delimiter=",");

	feat_embedded = TSNE(n_components=2).fit_transform(train_data)
	feat_embedded = np.asarray(feat_embedded)
	np.savetxt(train_tsne_2D_data, feat_embedded, delimiter=",");
