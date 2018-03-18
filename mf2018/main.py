# for generating and saving graphics without needing display
# 	NOTE: this MUST be at the beginning of the file
# 	see: https://stackoverflow.com/questions/15713279/calling-pylab-savefig-without-display-in-ipython
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from os import listdir # to get file names from a particular dir
import numpy as np
import classifier_helpers as ch
import matlab.engine

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
	tvt = ['train', 'valid', 'test'];

	# Read training and testing set
	for te_or_v_or_tr in tvt:
		# read boat set (TODO: print to make sure feature/labels appended correctly)
		filenames1 = read_filenames(datain_path+bh[B_LABEL]+'/'+te_or_v_or_tr+'/');
		labels1 = np.multiply(B_LABEL,np.ones(len(filenames1), dtype=int));
		# read heli set
		filenames2 = read_filenames(datain_path+bh[H_LABEL]+'/'+te_or_v_or_tr+'/');
		labels2 = np.multiply(H_LABEL,np.ones(len(filenames2), dtype=int))
		filenames = filenames1 + filenames2;
		labels = np.concatenate((labels1, labels2));
		# store to output dict
		dataset[te_or_v_or_tr+'_filenames'] = filenames;
		dataset[te_or_v_or_tr+'_labels'] = labels;
	
		print('Reading set: (' +te_or_v_or_tr+')');
		print('		'+str(len(filenames)) + ' files read');

	return dataset;

##############################
# MAIN CODE
##############################

### USER INPUT
# Define file use params
dataset_name = 'miniset'; # which dataset to use
navi_out = '../../'; # path to get out of github directory
vis_out = 'graphics/'; # path to output visuals

# Define audio params
sr = 44100; # sampling rate
clip_duration = 2; # audio clips will split into segments of this duration (in seconds)

# Define training params
max_epochs = 1; # maximum number of epochs for sgd

### START MATLAB ENGINE CODE
# This has to be in top-level file so it only has to start up once.. It takes a while
print('Starting MATLAB engine...');
mtl_eng = matlab.engine.start_matlab();
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
# TODO: see if this setup generates correctly, then add valid set before moving past file extraction
filepaths_dict = read_filenames_set(datain_path);

# ### TRAINING CLASSIFIER CODE
# Use training set to build classifier
train_filenames = filepaths_dict['train_filenames'];
train_labels = filepaths_dict['train_labels'];
# Use validation set to tune parameters
valid_filenames = filepaths_dict['valid_filenames'];
valid_labels = filepaths_dict['valid_labels'];
# Generate array containing classes
classes = np.array([B_LABEL, H_LABEL]);
print('*** Currently building classifier using training set');
clf_dict = ch.train_classifier(train_filenames, train_labels, valid_filenames, valid_labels, classes, max_epochs, sr, clip_duration, mtl_eng);
# Read output data and print
clf = clf_dict["clf"]; 
train_acc = clf_dict["train_acc"];
valid_acc = clf_dict["valid_acc"];
print(' ** final train accuracy: ');
print(train_acc);
print(' ** final validation accuracy: ');
print(valid_acc);

### EVALUATION CODE
# Use testing set to evaluate classifier
test_filenames = filepaths_dict['test_filenames'];
test_labels = filepaths_dict['test_labels'];
print('*** Currently evaluating classifier using testing set');
test_acc = ch.calc_accuracy(clf, test_filenames, test_labels, sr, clip_duration, mtl_eng);
print(' ** final test accuracy: ');
print(test_acc);


