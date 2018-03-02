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
# CHange so it doesn't chop, pass to feature gen instead
def import_data(filepath, vis_path, sr, clip_duration, plot=True):

	# Add subfolder to vis path
	vis_path = vis_path+'clips/';
	print(vis_path);

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

	# Loop and segment
	start_sec = 0;
	for n in range(num_clips):
		raw_audio = librosa_load_audio(filepath,sr,start_sec=start_sec,duration_sec=clip_duration);
		start_sec += clip_duration;

		# Generate magnitude spectrogram data
		specgram = mag_spec(raw_audio,sr);

		if plot:
			fig = plt.figure();
			librosa.display.waveplot(raw_audio, sr=sr);
			plt.title(filepath);
			plt.savefig(vis_path+filename+'_'+str(n)+'.png');
			plt.close(fig);
			fig = plt.figure();
			librosa.display.specshow(librosa.amplitude_to_db(specgram,ref=np.max),y_axis='log',x_axis='time');
			plt.title('Power spectrogram');
			plt.colorbar(format='%+2.0f dB');
			plt.tight_layout();
			plt.savefig(vis_path+filename+'_specgram_'+str(n)+'.png');
			plt.close(fig);

		

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
# mag_spec: generate magnitude spectrogram given audio
def mag_spec(raw_audio, sr):

	return librosa.stft(raw_audio);
	w_len = 0.02*sr;
	return librosa.stft(raw_audio,n_fft=w_len,hop_length=w_len/2,window='hamming');

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
	for filepath in filepaths_list:
		import_data(filepath, vis_path, sr, clip_duration, plot=True); # append a ,True if want to plot

