# gen_data.py: CODE PROTOTYPE - use main.py for the official version
# the data to use is defined in the main code segment
# this file reads in all the provided data and creates magnitude spectrogram features
# and outputs plots if requested

# for generating and saving graphics without needing display
# 	NOTE: this MUST be at the beginning of the file
# 	see: https://stackoverflow.com/questions/15713279/calling-pylab-savefig-without-display-in-ipython
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import librosa, librosa.display
import sys #system-specific parameters and functions
import numpy as np
import scipy as sp
import sklearn
from os import listdir # to get file names from a particular dir
import anurag_cnn.feat_extractor as cnn # to import cnn
from anurag_cnn.feat_extractor import usegpu

##############################
# GLOBALS
##############################
B_LABEL = 1;
H_LABEL = 2;

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
# EXTRACT AUDIO DATA SEGMENT
##############################
# import_data: top level call to import and segment files
# Inputs:
#	filepath - source of data
#	vis_path - path to output all graphs
#	sr - samp rate to convert all audio to
#	clip_duration - length in seconds to segment audio to
#	plot - toggle on visualization
# Output:
#
# Change so it doesn't chop, pass to feature gen instead
def import_data(key, filepaths_list, vis_path, sr, clip_duration, dataout_path, plot=False):

	rawaudio = 1;
	magspec = 1;

	for filepath in filepaths_list:

		# Output details
		filename = filepath.split('/');
		filename = filename[-1]; # get the last item on list, which should be filename
		print('Reading file: '+filename);

		# Check audio file's length in seconds (float) directly from file
		file_duration = librosa.get_duration(filename=filepath);
		print('	file duration: '+ str(file_duration));

		# Calculate number of clips audio will be segmented into
		num_clips = int(np.floor(file_duration/clip_duration));
		print('	num clips: '+ str(num_clips));

		# This section loops through the provided dataset purely for visualization of the raw data. 
		# Do not run if unneeded.
		if rawaudio:
			# Loop and segment
			start_sec = 0;
			for n in range(num_clips):
				raw_audio = librosa_load_audio(filepath,sr,start_sec=start_sec,duration_sec=clip_duration);
				start_sec += clip_duration;

				if plot:
					fig = plt.figure();
					librosa.display.waveplot(raw_audio, sr=sr);
					plt.title(filename);
					plt.savefig(vis_path+'rawaudio/'+filename+'_'+str(n)+'.png');
					plt.close(fig);

		# This section generates a magnitude spectrogram for each input audio clip.
		if magspec:
			magspec_out = run_mag_spec(filepath, sr, clip_duration, num_clips, vis_path, filename, plot);

			# Save with features along rows instead of columns
			magspec_out = np.transpose(magspec_out);
			n_features = np.shape(magspec_out)[0];

			fig = plt.figure();
			librosa.display.specshow(librosa.amplitude_to_db(magspec_out[300:2300,:], ref=np.max),y_axis='log',x_axis='time');
			plt.title(filename);
			plt.colorbar(format='%+2.0f dB');
			plt.tight_layout();
			plt.savefig(vis_path+'magspec/magspec_'+key+'.png');
			plt.close(fig);

			# Store outputs
			if ((key == 'boattest') or (key == 'boattrain')):
				# boat label: 0
				magspec_out = np.concatenate((np.zeros((n_features,1)),magspec_out), axis=1);
				np.savetxt(dataout_path+'magspec/magspec_'+key+'.csv', magspec_out, delimiter=',');
			if ((key == 'helitest') or (key == 'helitrain')):
				# heli label: 1
				magspec_out = np.concatenate((np.ones((n_features,1)),magspec_out), axis=1);
				np.savetxt(dataout_path+'magspec/magspec_'+key+'.csv', magspec_out, delimiter=',');


# librosa_load_audio: imports single audio file and return single raw data
def librosa_load_audio(filepath, sr, start_sec=0.0, duration_sec=None):
	# Import raw data
	try:
		raw_audio, sr = librosa.load(filepath, sr=sr, offset=start_sec, duration=duration_sec, mono=True);
	except:
		raise IOError('Give me an audio file which I can read!!');
	return raw_audio;

##############################
# FEATURE GENERATION SEGMENT
##############################
# run_mag_spec: loop through all files to call mag_spec
def run_mag_spec(filepath, sr, clip_duration, num_clips, vis_path, filename, plot=False):

	# Initialize matrix size
	magspec_shape = np.shape(mag_spec(librosa_load_audio(filepath,sr,start_sec=0,duration_sec=clip_duration), sr));
	magspec_out = np.empty((magspec_shape[0],0), dtype=float);

	# Initialize start time
	start_sec = 0;

	for n in range(num_clips):
		raw_audio = librosa_load_audio(filepath,sr,start_sec=start_sec,duration_sec=clip_duration);
		# raw_audio is np.array type

		# Generate magnitude spectrogram data
		specgram = mag_spec(raw_audio,sr);
		magspec_out = np.concatenate((magspec_out, specgram), axis=1);

		# Get start time for next clip
		start_sec += clip_duration;

		if plot:
			fig = plt.figure();
			librosa.display.specshow(librosa.amplitude_to_db(specgram, ref=np.max),y_axis='log',x_axis='time');
			plt.title(filename);
			plt.colorbar(format='%+2.0f dB');
			plt.tight_layout();
			plt.savefig(vis_path+'magspec/'+filename+'_'+str(n)+'.png');
			plt.close(fig);
	return magspec_out;

# gen_cnn: use a pre-trained CNN to extract features
def gen_cnn(raw_audio, use_gpu=True):
    
    features = cnn.main(raw_audio[0], use_gpu);
    
    return features;

# mag_spec: generate magnitude spectrogram given audio
def mag_spec(raw_audio, sr):

	return librosa.stft(raw_audio);
	w_len = 0.02*sr;
	return librosa.stft(raw_audio,n_fft=w_len,hop_length=round(w_len/2),window='hamming');

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
clip_duration = 1; # audio clips will split into segments of this duration (in seconds)

### EXTRACT FILENAMES CODE
# Construct pathways
datain_path = navi_out + 'data/' + dataset_name + '/'; # data input path
dataout_path = navi_out + 'data_out/' + dataset_name + '/'; # data exported here
vis_path = vis_out + dataset_name+'/'; # graphics exported here within repo

print('Reading data from: '+ datain_path);
print('Exporting data to: '+ dataout_path);

# Read filepaths set with keys 'boattest' 'boattrain' 'helitest' 'helitrain'
#	and entries as list of filepaths
filepaths_set = read_filenames_set(datain_path);
keys = ['boattest', 'boattrain', 'helitest', 'helitrain'];

### EXTRACT AUDIO CODE
for key in keys:
	filepaths_list = filepaths_set[key];
	print('*** Currently extracting for ' + key);
	import_data(key, filepaths_list, vis_path, sr, clip_duration, dataout_path, plot=True); # append a ,True if want to plot
	print('\n');
