import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from sklearn.manifold import TSNE
import numpy as np
import gmm_top_helpers as gmm_t
from sklearn.decomposition import PCA
from sklearn.externals import joblib
import svm_helper as svm_h

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
TRAIN_GMM = 0;
MAKE_SUPERVECTS = 0;  # run this every time new n-feats or ncl
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
feat_types = ['mod','spec','mel','spec-pca','mel-pca','mod-pca'];
feat_type = 'mel';

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
# Saving and reading trained PCA object
pcaobj_filepath = '/home/yichingl/mf/mf2018/mf_gmm/trained_svm_gmm/'+full_name+'.tsnepcaobj';


# Load supervects

## GMM
GMM = 1;
if GMM:
	# pack shared parameters
	data_params = {'n_feats': n_feats,'n_clusters': n_clusters, 'feat_type': feat_type};
	labels = {'BOAT_N': BOAT_N, 'HELI_N': HELI_N};
	# read test files
	path_params = {'boat_train_ctl_name': boat_train_ctl_name, 'heli_train_ctl_name': heli_train_ctl_name, 'ctlset_path': ctlset_path};
	train_data = gmm_t.read_train_data(data_params, path_params, labels);
	feat_read = train_data['train_feat_read'];
	labels_train = train_data['labels_train'];

## cnn 
else:
	# read control file for boat training data, create label vector, read group data
	filename = ctlset_path+boat_train_ctl_name;
	train_feat_read_boat = svm_h.read_cnn_features(filename+'.train.ctl',feat_type);
	labels_boat = svm_h.create_labels(train_feat_read_boat, BOAT_N);
	subgroups_boat = svm_h.read_subgroups(filename+'.train_subgroups.ctl');

	# read control file for heli training data, create label vector, read group data
	filename = ctlset_path+heli_train_ctl_name;
	train_feat_read_heli = svm_h.read_cnn_features(filename+'.train.ctl',feat_type);
	labels_heli = svm_h.create_labels(train_feat_read_heli, HELI_N);
	subgroups_heli = svm_h.read_subgroups(filename+'.train_subgroups.ctl');

	# concatenate boat and heli data
	feat_read = np.vstack((train_feat_read_boat,train_feat_read_heli));
	labels_train = np.concatenate((labels_boat,labels_heli));






# fit pca
pca = PCA(n_components = 50);
pca.fit(feat_read);
var = pca.explained_variance_ratio_;
var1 = np.cumsum(np.round(var, decimals=4)*100);
print(var1);

# 	# Save PCA transform obj
# 	joblib.dump(pca, pcaobj_filepath); 
# else:
# 	pca = joblib.load(pcaobj_filepath);

pca_feat = pca.transform(feat_read);
feat_embedded = TSNE(n_components=3).fit_transform(pca_feat)

if GMM:

	fig = plt.figure()
	ax = fig.add_subplot(111,projection='3d');
	plt.scatter(feat_embedded[:, 0], feat_embedded[:, 1],  feat_embedded[:, 2],c = labels_train);
	
	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')
	plt.savefig('../vis/tsne-gmm-spec-3d.png');
	pickle.dump(ax, file('tsne-gmm-spec-3d.pickle', 'w'))
	plt.close()

	plt.scatter(feat_embedded[:, 2], feat_embedded[:, 0], c = labels_train);
	plt.tight_layout();
	plt.savefig('../vis/tsne-gmm-spec-3d-zx.png');
	plt.close();

	plt.scatter(feat_embedded[:, 2], feat_embedded[:, 1], c = labels_train);
	plt.tight_layout();
	plt.savefig('../vis/tsne-gmm-spec-3d-zy.png');
	plt.close();


else:
	fig = plt.figure()
	ax = fig.add_subplot(111,projection='3d');
	plt.scatter(feat_embedded[:, 0], feat_embedded[:, 1],  feat_embedded[:, 2],c = labels_train);
	
	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')
	plt.savefig('../vis/tsne-cnn-spec-3d.png');
	pickle.dump(ax, file('tsne-cnn-spec-3d.pickle', 'w'))
	plt.close()

	plt.scatter(feat_embedded[:, 2], feat_embedded[:, 0], c = labels_train);
	plt.tight_layout();
	plt.savefig('../vis/tsne-cnn-spec-3d-zx.png');
	plt.close();

	plt.scatter(feat_embedded[:, 2], feat_embedded[:, 1], c = labels_train);
	plt.tight_layout();
	plt.savefig('../vis/tsne-cnn-spec-3d-zy.png');
	plt.close();



feat_embedded = TSNE(n_components=2).fit_transform(pca_feat)

if GMM:
	plt.scatter(feat_embedded[:, 0], feat_embedded[:, 1], c = labels_train);
	plt.tight_layout();
	plt.savefig('../vis/tsne-gmm-spec-2d.png');

else:
	plt.scatter(feat_embedded[:, 0], feat_embedded[:, 1], c = labels_train);
	plt.tight_layout();
	plt.savefig('../vis/tsne-cnn-spec-2d.png');



