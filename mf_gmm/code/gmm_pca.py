import numpy as np
from sklearn import preprocessing
import gmm_top_helpers as gmm_t
from sklearn.model_selection import GroupShuffleSplit
from sklearn import preprocessing
from sklearn.decomposition import PCA
import gmm_helper as gmm_h
from sklearn.externals import joblib

TEST = 0;

BOAT_N = 0;
HELI_N = 1;

def make_and_save_pcafeat(ctl_path, data_params, pca_obj, feat_type):

	# Unpack data parameters
	n_feats = data_params['n_feats'];
	n_clusters = data_params['n_clusters'];
	feat_type = data_params['feat_type'];

	with open(ctl_path) as f: # open boat,heli ctl file
		for line in f: # for each filepath, generate a supervect and save
			ctl_name = gmm_h.read_ctlname(line,feat_type);
			feat_read = gmm_h.read_specgram_feature(ctl_name);

			print('Reading ctl: ' + ctl_name);
			print('pre PCA transformed dimensions: ' + str(np.shape(feat_read)));
			pca_feat = pca_obj.transform(feat_read);
			print('PCA transformed dimensions: ' + str(np.shape(pca_feat)));

			# save means supervector
			pcafeat_out = ctl_name.replace(feat_type,feat_type+'-pca');
			with open(pcafeat_out, 'w') as f:
				np.savetxt(f,np.transpose(np.asarray(pca_feat)));
			print('Saving to: ' + pcafeat_out);


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
MAKE_SUPERVECTS = 0;  # run this every time new timeframe or ncl
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
n_clusters = 64;
max_iter = 5;
covariance_type = 'diag';

##############################
# SPECGRAM TYPE
##############################
feat_types = ['mod','spec','mel'];
feat_type = 'spec';

if feat_type == 'spec':
	#n_feats = 442;
	#n_timeframes = 200;

	#n_feats = 104;
	#n_timeframes = 849;

	n_feats = 133;
	n_timeframes = 669;

	pca_ncomp = 91; # cum var = 99.90999999999998

	#n_feats = 1025;
	#n_timeframes = 173;
	# if LEARN_PARAMS == 0:
	# 	best_kernel = 'linear';
	# 	best_C = 4096.0;
	# 	best_gamma: 'auto';

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
# Saving and reading trained PCA object
pcaobj_filepath = '/home/yichingl/mf/mf2018/mf_gmm/trained_svm_gmm/'+full_name+'.pcaobj';
##############################
# CNN PARAMETERS
##############################
# use xval to find best option out of the following
SVM_param_grid = [
{'C': np.logspace(-8.,12,num=11,base=2.0), 'gamma': ['auto'], 'kernel': ['linear']},
{'C': np.logspace(-8.,12,num=11,base=2.0), 'gamma': np.append(np.logspace(-5.,3,num=5,base=2.0),[1./n_feats]), 'kernel': ['rbf']},
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
# PCA
##############################
pca = PCA(n_components = pca_ncomp);

##############################
# MAIN CODE
##############################
# pack shared parameters
data_params = {'n_feats': n_feats,'n_clusters': n_clusters, 'feat_type': feat_type};
labels = {'BOAT_N': BOAT_N, 'HELI_N': HELI_N};

# load training data and fit pca
print('\n+++++++++++++++++++++++++++++++++++++++++++');
print('TRAINING GMM-UBM');
print('Using ctl names: ');
print(traingmm_ctl_names);
print('Savimg GMM-UBM to: ' + gmm_hdf5_filepath); 
path_params = {'traingmm_ctl_names': traingmm_ctl_names, 'gmm_hdf5_filepath': gmm_hdf5_filepath};
#gmm_t.train_gmm_mle(data_params, path_params, scaler);

# Unpack data parameters
n_feats = data_params['n_feats'];
n_clusters = data_params['n_clusters'];
feat_type = data_params['feat_type'];
# Unpack path parameters
traingmm_ctl_names = path_params['traingmm_ctl_names'];
gmm_hdf5_filepath = path_params['gmm_hdf5_filepath'];

# Initialize matrix for features
feat_all = np.zeros((0,n_feats));

# Read all ctl files
for ctl_group_name in traingmm_ctl_names:

	# Generate input ctl read path
	ctl_path = '../ctl_files/' + ctl_group_name +'.train.ctl';

	# Read ctl file
	read_feat = gmm_h.read_specgram_features(ctl_path, n_feats, feat_type);

	# Add to full matrix
	feat_all = np.vstack((feat_all, read_feat));

# Normalize data along feature axis
scaler.fit(feat_all);
feat_all = scaler.transform(feat_all);

# Perform PCA
pca.fit(feat_all);

# Load and transform training set
# Process training set (boat, heli)
for ctl_group_name_supervect in gensupervect_train_ctl_names: # boat, heli
	ctl_path = '../ctl_files/' + ctl_group_name_supervect +'.train.ctl';
	make_and_save_pcafeat(ctl_path, data_params, pca,feat_type);

# Process testing set (boat, heli)
for ctl_group_name_supervect in gensupervect_test_ctl_names: # boat, heli
	ctl_path = '../ctl_files/' + ctl_group_name_supervect +'.test.ctl';
	make_and_save_pcafeat(ctl_path, data_params, pca,feat_type);

# Save PCA transform obj
joblib.dump(pca, pcaobj_filepath); 
