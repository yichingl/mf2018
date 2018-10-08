import numpy as np
import gmm_top_helpers as gmm_t
import pandas as pd

TEST = 0;

BOAT_N = 0;
HELI_N = 1;

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

if TEST:
	print('TESTING');
	group_name = 'test';
	# For SVM
	boat_train_ctl_name = 'test.boat';
	heli_train_ctl_name = 'test.heli';
	boat_test_ctl_name = 'minitest.boatA';
	heli_test_ctl_name = 'minitest.heliA';
	# For GMM
	traingmm_ctl_names = ['test.boat','test.heli'];
	gensupervect_train_ctl_names = [boat_train_ctl_name, heli_train_ctl_name];
	gensupervect_test_ctl_names = [boat_test_ctl_name, heli_test_ctl_name];

##############################
# TOGGLE TASKS
##############################
TRAIN_GMM = 1;
MAKE_SUPERVECTS = 1;  # run this every time new n-feats or ncl
LEARN_PARAMS = 1; 
TRAIN_SVM = 1; # must use learn_params together? for now
CLASSIFY = 1;

##############################
# AUDIO PARAMETERS
##############################
sr = 44100; # sampling rate
clip_duration = 2; # audio clips will split into segments of this duration (in seconds)

##############################
# GMM PARAMETERS
##############################
n_clusters = 16;
max_iter = 5;
covariance_type = 'diag';

##############################
# SPECGRAM TYPE
##############################
feat_types = ['mod','spec','mel','spec-pca','mel-pca','mod-pca'];
feat_type = 'spec';

if feat_type == 'spec':
	#n_feats = 442;
	#n_timeframes = 200;

	#n_feats = 104;
	#n_timeframes = 849;

	n_feats = 133;
	n_timeframes = 669;

	#n_feats = 1025;
	#n_timeframes = 173;
	# if LEARN_PARAMS == 0:
	# 	best_kernel = 'linear';
	# 	best_C = 4096.0;
	# 	best_gamma: 'auto';
elif feat_type == 'spec-pca':
	n_feats = 91;
	n_timeframes = 669;

elif feat_type == 'mel':
	n_feats = 128;
	n_timeframes = 173;

	# if LEARN_PARAMS == 0:
	# 	best_kernel = 'rbf';
	# 	best_C = 4096.0;
	# 	best_gamma = 0.001;

	# n_feats = 133;
	# n_timeframes = 669;
elif feat_type == 'mod':
	n_feats = 128;
	n_timeframes = 320;
	# if LEARN_PARAMS == 0:
	# 	best_kernel = 'rbf';
	# 	best_C = 1024.0;
	# 	best_gamma = 'auto';

##############################
# FULL_NAME
##############################
full_name = group_name+'.'+feat_type+'.'+str(n_clusters)+'n_cl'+'.'+str(n_feats)+'n_feats';

##############################
# I/O PATHS
##############################
# Location of ctl files
ctlset_path = '../ctl_files/'
# Reading spectrogram path
feat_out = '/home/yichingl/mf/feat_out/';
# Writing supervect path
supervect_out = '/home/yichingl/mf/supervect_out/';
# Saving and reading GMM-UBM path
gmm_hdf5_filepath = '/home/yichingl/mf/mf2018/mf_gmm/trained_gmms/'+full_name+'.GMM_UBM.hdf5';
# Saving and reading trained SVM
trained_svm_filepath = '/home/yichingl/mf/mf2018/mf_gmm/trained_svm_gmm/'+full_name+'.GMM_SVM.sav';
# Saving and reading trained SVM scaler object
trained_svm_scalerobj_filepath = '/home/yichingl/mf/mf2018/mf_gmm/trained_svm_gmm/'+full_name+'.GMM_SVM.scalerobj';
# Pretty print results filepath
pretty_results_filepath = '../new_gmm_results/' + full_name +'.results.txt';
# Best params
best_params_filepath ='../new_gmm_results/' + full_name +'.best_params.txt';
# Dataframe output path
train_dataframe_out_filepath = supervect_out+full_name+'.train_dataframe.csv';
test_dataframe_out_filepath = supervect_out+full_name+'.test_dataframe.csv';

# pack shared parameters
data_params = {'n_feats': n_feats,'n_clusters': n_clusters, 'feat_type': feat_type};
labels = {'BOAT_N': BOAT_N, 'HELI_N': HELI_N};

# convert individual specgram files into a single pandas dataframe
print('\n+++++++++++++++++++++++++++++++++++++++++++');
print('LOADING TRAINING DATA FOR SVM');
print('Using boat train ctl name: ');
print(boat_train_ctl_name);
print('Using heli train ctl name: ');
print(heli_train_ctl_name);

# read test files
path_params = {'boat_train_ctl_name': boat_train_ctl_name, 'heli_train_ctl_name': heli_train_ctl_name, 'ctlset_path': ctlset_path};
train_data = gmm_t.read_train_data(data_params, path_params, labels);

# extract data from dict
train_feat_read = train_data['train_feat_read'];
print('Train Data: ');
print(train_feat_read);
labels_train = train_data['labels_train'];
subgroups = train_data['subgroups'];
# generate dataframe row and col labels
sample_inds = np.arange(np.shape(train_feat_read)[0]);
feat_names = ['chann_'+str(n) for n in np.arange(np.shape(train_feat_read)[1])];
# create dataframe object
train_dataframe = pd.DataFrame(data=train_feat_read,index=sample_inds, columns=feat_names);
print('Train DataFrame after initial setup:');
print(train_dataframe);
# add data lables and subgroups
train_dataframe = train_dataframe.assign(labels_train=labels_train,subgroups=subgroups);
print('Train DataFrame after assigning labels and subgroups:');
print(train_dataframe);
# save pandas array
train_dataframe.to_csv(train_dataframe_out_filepath);


# read test files
path_params = {'boat_test_ctl_name': boat_test_ctl_name, 'heli_test_ctl_name': heli_test_ctl_name, 'ctlset_path': ctlset_path};
test_data = gmm_t.read_test_data(data_params, path_params, labels);

# extract data from dict
test_feat_read = test_data['test_feat_read'];
print('Test Data: ');
print(test_feat_read);
labels_test = test_data['labels_test'];
# generate dataframe row and col labels
sample_inds = np.arange(np.shape(test_feat_read)[0]);
feat_names = ['chann_'+str(n) for n in np.arange(np.shape(test_feat_read)[1])];
# create dataframe object
test_dataframe = pd.DataFrame(data=test_feat_read,index=sample_inds, columns=feat_names);
print('Test DataFrame after initial setup:');
print(train_dataframe);
# add data lables and subgroups
test_dataframe = test_dataframe.assign(labels_test=labels_test);
print('Test DataFrame after assigning labels:');
print(train_dataframe);
# save pandas array
test_dataframe.to_csv(test_dataframe_out_filepath);