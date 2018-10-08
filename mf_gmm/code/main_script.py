import numpy as np
from sklearn import preprocessing
import gmm_top_helpers as gmm_t
from sklearn.model_selection import GroupShuffleSplit
from sklearn import preprocessing

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
LEARN_PARAMS = 0; 
TRAIN_SVM = 0; # must use learn_params together? for now
CLASSIFY = 0;

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
feat_type = 'mod';

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
gmm_hdf5_filepath = '../../../trained_gmms/'+full_name+'.GMM_UBM.hdf5';
# Saving and reading trained SVM
trained_svm_filepath = '/home/yichingl/mf/mf2018/mf_gmm/trained_svm_gmm/'+full_name+'.GMM_SVM.sav';
# Saving and reading trained SVM scaler object
trained_svm_scalerobj_filepath = '/home/yichingl/mf/mf2018/mf_gmm/trained_svm_gmm/'+full_name+'.GMM_SVM.scalerobj';
# Pretty print results filepath
pretty_results_filepath = '../new_gmm_results/' + full_name +'.results.txt';
# Best params
best_params_filepath ='../new_gmm_results/' + full_name +'.best_params.txt';
##############################
# CNN PARAMETERS
##############################
# use xval to find best option out of the following
SVM_param_grid = [
{'C': np.logspace(-8.,12,num=11,base=2.0), 'gamma': ['auto'], 'kernel': ['linear']},
{'C': np.logspace(-8.,12,num=11,base=2.0), 'gamma': np.append(np.logspace(-8.,3,num=7,base=2.0),[1./n_feats]), 'kernel': ['rbf']},
]
xval_test_size = 1./3;
xval_num_splits = 3;
split_obj = GroupShuffleSplit(test_size=xval_test_size, n_splits=xval_num_splits);

##############################
# NORMALIZATION METHOD
##############################
norm_method = 'minmax'; # 'standard'
if norm_method == 'minmax':
	scaler = preprocessing.MinMaxScaler();
elif norm_method == 'standard':
	scaler = preprocessing.StandardScaler();

##############################
# MAIN CODE
##############################
# pack shared parameters
data_params = {'n_feats': n_feats,'n_clusters': n_clusters, 'feat_type': feat_type};
labels = {'BOAT_N': BOAT_N, 'HELI_N': HELI_N};

# train and save GMM-UBM
if TRAIN_GMM:
	print('\n+++++++++++++++++++++++++++++++++++++++++++');
	print('TRAINING GMM-UBM');
	print('Using ctl names: ');
	print(traingmm_ctl_names);
	print('Savimg GMM-UBM to: ' + gmm_hdf5_filepath); 
	path_params = {'traingmm_ctl_names': traingmm_ctl_names, 'gmm_hdf5_filepath': gmm_hdf5_filepath};
	gmm_t.train_gmm_mle(data_params, path_params, scaler);

# load GMM-UBM and MAP adapt all train and test files that will be used by SVM
if MAKE_SUPERVECTS:
	print('\n+++++++++++++++++++++++++++++++++++++++++++');
	print('MAKE SUPERVECTS');
	print('Using train ctl names: ');
	print(gensupervect_train_ctl_names);
	print('Using test ctl names: ');
	print(gensupervect_test_ctl_names);
	print('Loading GMM-UBM from: ' + gmm_hdf5_filepath);
	path_params = {'gensupervect_train_ctl_names': gensupervect_train_ctl_names,'gensupervect_test_ctl_names': gensupervect_test_ctl_names, 'gmm_hdf5_filepath': gmm_hdf5_filepath};
	gmm_t.make_supervects(data_params, path_params);

# read training data
if LEARN_PARAMS or TRAIN_SVM:
	print('\n+++++++++++++++++++++++++++++++++++++++++++');
	print('LOADING TRAINING DATA FOR SVM');
	print('Using boat train ctl name: ');
	print(boat_train_ctl_name);
	print('Using heli train ctl name: ');
	print(heli_train_ctl_name);
	print('Loading GMM-UBM from: ' + gmm_hdf5_filepath);
	path_params = {'boat_train_ctl_name': boat_train_ctl_name, 'heli_train_ctl_name': heli_train_ctl_name, 'ctlset_path': ctlset_path};
	train_data = gmm_t.read_train_data(data_params, path_params, labels);


# learn best parameters for SVM, then save
if LEARN_PARAMS:
	print('\n+++++++++++++++++++++++++++++++++++++++++++');
	print('LEARNING SVM PARAMS');
	print('Xval test size: '+ str(xval_test_size));
	print('Xval num splits: ' + str(xval_num_splits));
	print('Saving pretty print to: ' + pretty_results_filepath);
	print('Saving best params to: '+ best_params_filepath); # RESUME HERE: make this global
	xval_params = {'xval_test_size': xval_test_size, 'xval_num_splits': xval_num_splits};
	gmm_t.xval_params(train_data, SVM_param_grid, xval_params, full_name,pretty_results_filepath,best_params_filepath,split_obj,scaler);

# use best parameter to train and save an SVM
if TRAIN_SVM:

	path_params = {'full_name': full_name, 'trained_svm_filepath': trained_svm_filepath, 'trained_svm_scalerobj_filepath': trained_svm_scalerobj_filepath};
	#SVM_params = {'best_kernel': best_kernel, 'best_C': best_C, 'best_gamma': best_gamma};
	SVM_params = gmm_t.read_SVM_param(best_params_filepath,n_feats);

	gmm_t.train_svm(train_data, path_params, SVM_params, scaler);

# load test data and classify
if CLASSIFY:
	# read test files
	path_params = {'boat_test_ctl_name': boat_test_ctl_name, 'heli_test_ctl_name': heli_test_ctl_name, 'ctlset_path': ctlset_path};
	test_data = gmm_t.read_test_data(data_params, path_params, labels);

	# classify
	path_params = {'full_name': full_name, 'trained_svm_filepath': trained_svm_filepath, 'trained_svm_scalerobj_filepath': trained_svm_scalerobj_filepath};
	gmm_t.classify_data(test_data, path_params, scaler);