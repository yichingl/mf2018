import numpy as np;
import sklearn.metrics as sklmetrics

##############################
# READ DATA
##############################
# read control file and return data matrix
def read_cnn_features(filename, feat_type):
	train_feat_read = np.zeros((0,1024));
	with open(filename) as f:
		# extract file name, and start and end vectors
		for line in f:
			ctl_name = line.rstrip('\n').replace('feat_out','CNN_feat_out');
			ctl_name = ctl_name.replace('spec',feat_type);

			# read data from control file into 2D matrix form
			feat_read = np.loadtxt(ctl_name);

			# concatenate
			train_feat_read = np.vstack((train_feat_read,feat_read));
	return np.absolute(train_feat_read);

# read control file and return data matrix
def read_gmm_features(filename, feat_type, feat_size):
	train_feat_read = np.zeros((0,feat_size));
	with open(filename) as f:
		# extract file name, and start and end vectors
		for line in f:
			ctl_name = line.rstrip('\n').replace('feat_out','supervect_out');
			ctl_name = ctl_name.replace('spec',feat_type);

			# read data from control file into 2D matrix form
			feat_read = np.loadtxt(ctl_name);

			# concatenate
			train_feat_read = np.vstack((train_feat_read,feat_read)); 

	return np.absolute(train_feat_read);

# read data matrix size and generate corresponding label vector
def create_labels(data_mat, label):
	labels = np.ones(np.shape(data_mat)[0],dtype=int)*label;
	labels = labels.astype(int);
	return labels;

# read subgroup data from subgroup filename
def read_subgroups(filename):
	with open(filename) as f:
		subgroups = f.read().splitlines();
	return subgroups;

##############################
# TRAIN SVM
##############################
def calc_accuracy(true_labels,predicted_labels):
	return sklmetrics.accuracy_score(true_labels,predicted_labels);

def calc_f1(true_labels,predicted_labels):
	return sklmetrics.f1_score(true_labels,predicted_labels);

def calc_mcc(true_labels,predicted_labels):
	return sklmetrics.matthews_corrcoef(true_labels,predicted_labels);

def calc_precision(true_labels,predicted_labels):
	return sklmetrics.precision_score(true_labels,predicted_labels);

def calc_recall(true_labels,predicted_labels):
	return sklmetrics.recall_score(true_labels,predicted_labels);

def calc_cfmat(true_labels,predicted_labels):
	return sklmetrics.confusion_matrix(true_labels,predicted_labels);

