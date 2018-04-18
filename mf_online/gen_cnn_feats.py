

# pregenerate and store cnn features

# for generating and saving graphics without needing display
# 	NOTE: this MUST be at the beginning of the file
# 	see: https://stackoverflow.com/questions/15713279/calling-pylab-savefig-without-display-in-ipython

from os import listdir # to get file names from a particular dir
import numpy as np
import features_helpers as fh

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

# read_featnames: read all feat names from given folder and return list
def read_featnames(inpath, outpath):
	filenames = [outpath + f.replace(".wav", ".csv") for f in listdir(inpath) if (f[-4:] == ".wav")]; 
	return filenames;

# read_filenames_set: read in all files under train,test types in a dict with 
#	keys 'train' 'test'
#	and entries as lists of filenames
def read_filenames_set(datain_path,featpath):
	dataset = {};
	bh = ['boat', 'heli'];
	tvt = ['train', 'valid', 'test'];

	# Read training and testing set
	for te_or_v_or_tr in tvt:
		# read boat set (TODO: print to make sure feature/labels appended correctly)
		filenames1 = read_filenames(datain_path+bh[B_LABEL]+'/'+te_or_v_or_tr+'/');
		featnames1 = read_featnames(datain_path+bh[B_LABEL]+'/'+te_or_v_or_tr+'/', featpath+bh[B_LABEL]+'/'+te_or_v_or_tr+'/');
		# read heli set
		filenames2 = read_filenames(datain_path+bh[H_LABEL]+'/'+te_or_v_or_tr+'/');
		featnames2 = read_featnames(datain_path+bh[H_LABEL]+'/'+te_or_v_or_tr+'/', featpath+bh[H_LABEL]+'/'+te_or_v_or_tr+'/');
		filenames = filenames1 + filenames2;
		featnames = featnames1 + featnames2;
		# store to output dict
		dataset[te_or_v_or_tr+'_filenames'] = filenames;
		dataset['feat_'+te_or_v_or_tr+'_filenames'] = featnames;
	
		print('Reading set: (' +te_or_v_or_tr+')');
		print('		'+str(len(filenames)) + ' files read');


	return dataset;

##############################
# MAIN CODE
##############################

### USER INPUT
# Define file use params
dataset_name = 'audio_used'; # which dataset to use
navi_out = '../../'; # path to get out of github directory
vis_out = 'graphics/'; # path to output visuals'

# Define audio params
sr = 44100; # sampling rate
clip_duration = 2; # audio clips will split into segments of this duration (in seconds)

# Define training params
max_epochs = 11; # maximum number of epochs for sgd
min_epochs = 5; # minimum number of epochs for sgd

### START MATLAB ENGINE CODE

### EXTRACT FILENAMES CODE
# Construct pathways
datain_path = navi_out + 'data/' + dataset_name + '/'; # data input path
dataout_path = navi_out + 'data_out/' + dataset_name + '/'; # data exported here
featpath = navi_out + 'data/' + dataset_name + '/feat_out/'; # features exported here
vis_path = vis_out + dataset_name+'/'; # graphics exported here within repo

print('Reading data from: '+ datain_path);
print('Exporting data to: '+ dataout_path);

# Read filepaths_dict dict with keys 'train_filenames' 'train_labels' 'valid_filenames' 'valid_labels' 'test_filenames' 'test_labels'
#	and entries as list of filepaths (train/test_filenames)
#		and mats of labels (train/test_labels)
# TODO: see if this setup generates correctly, then add valid set before moving past file extraction
filepaths_dict = read_filenames_set(datain_path,featpath);

# ### TRAINING CLASSIFIER CODE
# Use training set to build classifier
train_filenames = filepaths_dict['train_filenames'];
train_featnames = filepaths_dict['feat_train_filenames'];
fh.store_cnn(train_filenames, train_featnames, sr, clip_duration);
# Use validation set to tune parameters
valid_filenames = filepaths_dict['valid_filenames'];
valid_featnames = filepaths_dict['feat_valid_filenames'];
fh.store_cnn(valid_filenames, valid_featnames, sr, clip_duration);
# Use test set to tune parameters
test_filenames = filepaths_dict['test_filenames'];
test_featnames = filepaths_dict['feat_test_filenames'];
fh.store_cnn(test_filenames, test_featnames, sr, clip_duration);