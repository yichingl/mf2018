# for generating and saving graphics without needing display
# 	NOTE: this MUST be at the beginning of the file
# 	see: https://stackoverflow.com/questions/15713279/calling-pylab-savefig-without-display-in-ipython
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import librosa, librosa.display
import sys #system-specific parameters and functions
import numpy as np
from os import listdir # to get file names from a particular dir


##############################
# IMPORT DATAPATHS SEGMENT
##############################

# read_filenames: read all file names from given folder and return list
def read_filenames(path):
	filenames = [path + f for f in listdir(path) if (f[-4:] == ".wav")];  # include path
	#filenames = [f for f in listdir(path) if (f[-4:] == ".wav")]; # don't include path
	return filenames;

# read_filenames_set: read in all files under their respective types in a dict with 
#	keys 'boattest' 'boattrain' 'helitest' 'helitrain'
#	and entries as lists of filenames
def read_filenames_set(datain_path):
	dataset = {};
	bh = ['boat', 'heli'];
	tt = ['test', 'train'];
	for boat_or_heli in bh:
		for test_or_train in tt:
			filenames = read_filenames(datain_path+boat_or_heli+'/'+test_or_train+'/');
			dataset[boat_or_heli+test_or_train] = filenames;
			print('Reading set: (' +boat_or_heli+', '+test_or_train+')');
			print('		'+str(len(filenames)) + ' files read');
	return dataset;

##############################
# EXTRACT AUDIO SEGMENT
##############################
# import_wav: imports audio file as raw data
def import_wav(filepath, vis_path, sr, plot=False):
    # Import raw data
    try:
        raw_data, sr = librosa.load(filepath, sr=sr, mono=True);
    except:
        raise IOError('Give me an audio file which I can read!!');
    
    # Only use one channel
    #if len(raw_data[0].shape) > 1:
    #    raw_data = (raw_data[0][0], raw_data[1])
    
    if plot:
    	filename = filepath.split('/');
    	filename = filename[-1]; # get the last item on list, which should be filename
    	print('reading and saving figure for: '+filename);
        plt.figure();
        librosa.display.waveplot(raw_data, sr=sr);
        plt.title(filepath);
        plt.savefig(vis_path+filename+'read_audio.png');
        #plt.show()

    return raw_data
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
clip_duration = None # audio clips will split into segments of this duration

### EXTRACT FILENAMES CODE
# Construct pathways
datain_path = navi_out + 'data/' + dataset_name + '/'; # data input path
dataout_path = navi_out + 'data_out/' + dataset_name + '/'; # data exported here
vis_path = vis_out + dataset_name; # graphics exported here within repo

print('Reading data from: '+ datain_path);
print('Exporting data to: '+ dataout_path);

# Read filepaths set with keys 'boattest' 'boattrain' 'helitest' 'helitrain'
#	and entries as list of filepaths
filepaths_set = read_filenames_set(datain_path);
keys = ['boattest', 'boattrain', 'helitest', 'helitrain'];

### EXTRACT AUDIO CODE


# y,sr = librosa.load(librosa.util.example_audio_file(),sr=sr,mono=True);
# plt.figure()
# librosa.display.waveplot(y, sr=sr)
# plt.title('Monophonic')

for key in keys:
	filepaths_list = filepaths_set[key];
	for filepath in filepaths_list:
		import_wav(filepath, vis_path, sr, 1);

