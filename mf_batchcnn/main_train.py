# training main: Generate datapaths, then fine tune transCNN for each specgram type.  
#	Pass train data through each CNN and extract features.  Perhaps see the accuracy using individual CNNs.
#	Then train an SVM using these features.  Store each model. 

# for generating and saving graphics without needing display
# 	NOTE: this MUST be at the beginning of the file
# 	see: https://stackoverflow.com/questions/15713279/calling-pylab-savefig-without-display-in-ipython
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

from os import listdir # to get file names from a particular dir
import specCNN_helpers as sCNNh
import matlab.engine
import numpy as np
import random # for shuffling

##############################
# GLOBALS
##############################
B_LABEL = 0;
H_LABEL = 1;

##############################
# IMPORT DATAPATHS SEGMENT
##############################

# read_filenames: read all file names from given folder and return list
def read_filenames(path):
	filenames = [path + f for f in listdir(path) if (f[-4:] == ".wav")];  # include path
	#filenames = [f for f in listdir(path) if (f[-4:] == ".wav")]; # don't include path
	return filenames;

# read_filenames_set: read in all files under train,test types in a dict with 
#	keys 'train' 'test'
#	and entries as lists of filenames
def read_filenames_set(datain_path):
	dataset = {};
	bh = ['boat', 'heli'];
	tt = ['train'];

	# Read training and testing set
	for tr_or_te in tt:
		# read boat set (TODO: print to make sure feature/labels appended correctly)
		filenames1 = read_filenames(datain_path+bh[B_LABEL]+'/'+tr_or_te+'/');
		labels1 = np.multiply(B_LABEL,np.ones(len(filenames1), dtype=int));
		# read heli set
		filenames2 = read_filenames(datain_path+bh[H_LABEL]+'/'+tr_or_te+'/');
		labels2 = np.multiply(H_LABEL,np.ones(len(filenames2), dtype=int))
		filenames = filenames1 + filenames2;
		labels = np.concatenate((labels1, labels2));
		# store to output dict
		dataset[tr_or_te+'_filenames'] = filenames;
		dataset[tr_or_te+'_labels'] = labels;
	
		print('Reading set: (' +tr_or_te+')');
		print('		'+str(len(filenames)) + ' files read');

	return dataset;

# data_class holds train filenames, labels
class DataClass:

	def __init__(self, filenames, labels, batch_size):

		# filename and label data
		self.filenames = filenames;
		self.labels = labels;
		self.n_examples = len(filenames);

		# define how many samples to use each time
		self.batch_size = batch_size;

		# initialize n_shuffle
		self.n_shuffle = np.zeros(self.n_examples);

		# initialize batch_done bool that alerts us that there are no more examples to use
		self.batch_done = 0;
		self.cur_batch_ind = 0;

		# initialize cur batch dataset
		self.cur_batch_filenames = 0;
		self.cur_batch_labels = 0;

	def reset_shuffle(self):
		# shuffle data
		self.n_shuffle = np.arange(self.n_examples,dtype=np.int32);
		np.random.shuffle(self.n_shuffle);

		# reset batch_done and start index
		self.batch_done = 0;
		self.cur_batch_ind = np.int32(0);

	def make_batch(self):
		# load variables for convenience
		batch_size = self.batch_size;
		cur_batch_ind_start = self.cur_batch_ind;
		cur_batch_ind_end = cur_batch_ind_start+batch_size;
		print('current batch ind start: ' + str(cur_batch_ind_start));
		print('current batch size: ' + str(batch_size));
		print('current batch ind end: ' + str(cur_batch_ind_end));
		n_examples = self.n_examples;

		# if at the end of list, reset
		if cur_batch_ind_end >= n_examples:
			cur_batch_ind_end = n_examples;
			self.batch_done = 1;

		cur_batch_inds = self.n_shuffle[cur_batch_ind_start:cur_batch_ind_end];
		cur_batch_filenames = [];
		for cur_batch_ind in cur_batch_inds:
			cur_batch_filenames.append(self.filenames[cur_batch_ind]);
		self.cur_batch_labels = self.labels[cur_batch_inds];
		self.cur_batch_filenames = cur_batch_filenames;
		self.cur_batch_ind = cur_batch_ind_end;

##############################
# MAIN CODE: Training model
##############################

### USER INPUT
# Define file use params
dataset_name = 'mini'; # which dataset to use
navi_out = '../../'; # path to get out of github directory
vis_out = 'graphics/'; # path to output visuals

# Batch size (number of files in a batch)
batch_size = 5;

# Define audio params
sr = 44100; # sampling rate
clip_duration = 2; # audio clips will split into segments of this duration (in seconds)

# Define training params
cnn_max_epochs = 1;
svm_max_epochs = 1;

# Which specgrams?
use_magspec = 1;
train_magspec = 1;
use_modspec = 0;
train_modspec = 1;
use_melspec = 0;
train_melspec = 1;

### START MATLAB ENGINE CODE
# This has to be in top-level file so it only has to start up once.. It takes a while
print('Starting MATLAB engine...');
#mtl_eng = matlab.engine.start_matlab();
mtl_eng = matlab.engine.start_matlab('-nojvm');
mtl_eng.addpath(mtl_eng.genpath("matlab_toolboxes/"));
print('Done importing MATLAB engine...');

### EXTRACT FILENAMES CODE
# Construct pathways
datain_path = navi_out + 'data/' + dataset_name + '/'; # data input path
dataout_path = navi_out + 'data_out/' + dataset_name + '/'; # data exported here
vis_path = vis_out + dataset_name+'/'; # graphics exported here within repo

print('Reading data from: '+ datain_path);
print('Exporting data to: '+ dataout_path);

# Read filepaths_dict dict with keys 'train_filenames' 'train_labels' 'valid_filenames' 'valid_labels' 'test_filenames' 'test_labels'
#	and entries as list of filepaths (train/test_filenames)
#		and mats of labels (train/test_labels)
filepaths_dict = read_filenames_set(datain_path);

### TRAINING CLASSIFIER CODE
# Generate training set filepaths to be used to build classifier
train_filenames = filepaths_dict['train_filenames'];
train_labels = filepaths_dict['train_labels'];
data_class = DataClass(train_filenames,train_labels,batch_size);
print(data_class.filenames);
data_class.reset_shuffle();
print(data_class.n_shuffle);
while(not data_class.batch_done):
	data_class.make_batch();
	print(data_class.cur_batch_filenames);
	print(data_class.cur_batch_labels);
	print(data_class.cur_batch_ind);
data_class.reset_shuffle();
print(data_class.n_shuffle);
while(not data_class.batch_done):
	data_class.make_batch();
	print(data_class.cur_batch_filenames);
	print(data_class.cur_batch_labels);
	print(data_class.cur_batch_ind);

# Fine tune transCNN for each spectrogram type
classes = np.array([B_LABEL, H_LABEL]);
if use_magspec:
	if train_magspec:
		print('*** Currently training magnitude spectrogram CNN classifier using training set');
		magspec_CNN_clf = ch.finetuneCNN(train_filenames, train_labels, classes, cnn_max_epochs, sr, clip_duration, mtl_eng, 'magspec', batch_size);
	#else:
		# load the model
	# Pass training data back through CNN clfs and extract features
	magspec_CNN_feats = ch.extract_feats(magspec_CNN_feats,train_filenames);

if use_modspec:
	if train_modspec:
		print('*** Currently training modulation spectrogram CNN classifier using training set');
		modspec_CNN_clf = ch.finetuneCNN(train_filenames, train_labels, classes, cnn_max_epochs, sr, clip_duration, mtl_eng, 'modspec', batch_size);
	#else:
		# load the model
	# Pass training data back through CNN clfs and extract features
	modspec_CNN_feats = ch.extract_feats(modspec_CNN_feats,train_filenames);
if use_modspec:
	if train_modspec:
		print('*** Currently training mel spectrogram CNN classifier using training set');
		melspec_CNN_clf = ch.finetuneCNN(train_filenames, train_labels, classes, cnn_max_epochs, sr, clip_duration, mtl_eng, 'melspec', batch_size);
	#else:
		# load the model
	# Pass training data back through CNN clfs and extract features
	melspec_CNN_feats = ch.extract_feats(melspec_CNN_feats,train_filenames);
# Save each fine-tuned model

# Concat features
#CNN_feats = magspec_CNN_feats;
#CNN_feats = np.hstack(magspec_CNN_feats,modspec_CNN_feats,melspec_CNN_feats);

# Train RBF SVM on CNN feats
#SVM_clf = ch.trainSVM(CNN_feats, labels);
# Save SVM model






# Get training error? on subset?
