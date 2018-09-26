
from os import listdir # to get file names from a particular dir
import numpy as np
import librosa
import matlab.engine


##############################
# IMPORT DATAPATHS SEGMENT
##############################

# read_filenames: read all file names from given folder and return list
def read_filenames(path):
	filenames = [path + f for f in listdir(path) if (f[-4:] == ".wav")];  # include path
	#filenames = [f for f in listdir(path) if (f[-4:] == ".wav")]; # don't include path
	return filenames;

# make_out_filenames: read all file names from given folder and return list of output filenames
def make_out_filenames(in_path,out_path):
	#filenames = [path + f for f in listdir(path) if (f[-4:] == ".wav")];  # include path
	filenames = [out_path + f[0:-4] for f in listdir(in_path) if (f[-4:] == ".wav")]; # don't include path
	return filenames;

# read_filenames_set: read in all files under train,test types in a dict with 
#	keys 'train' 'test'
#	and entries as lists of filenames
def read_filenames_set(datain_path, dataout_path):
	dataset = {};

	# read set (TODO: print to make sure feature appended correctly)
	in_filenames = read_filenames(datain_path);
	out_filenames = make_out_filenames(datain_path,dataout_path);
		
	# store to output dict
	dataset['in_filenames'] = in_filenames;
	dataset['out_filenames'] = out_filenames;
	
	return dataset;

##############################
# IMPORT AUDIO
##############################
# load_audio: imports single audio file and return single raw data
def load_audio(filepath, sr, start_sec=0.0, duration_sec=None):
	# Import raw data
	try:
		raw_audio, sr = librosa.load(filepath, sr=sr, offset=start_sec, duration=duration_sec, mono=True);
	except:
		raise IOError('Give me an audio file which I can read!!');
	return raw_audio;

##############################
# OBTAIN DIRECTORIES USED
##############################
def read_train_file(directory_list_filepath):
	COPY = False;
	with open(directory_list_filepath, 'r') as f:
		for line in f:
			if line.strip() == "TRAIN_DIR":
				TRAIN_DIRS = [];
				COPY = True;
			elif line.strip() == "":
				COPY = False;
			elif COPY:
				TRAIN_DIRS.append(line.strip());
	return TRAIN_DIRS;
def read_test_file(directory_list_filepath):
	COPY = False;
	with open(directory_list_filepath, 'r') as f:
		for line in f:
			if line.strip() == "TEST_DIR":
				TEST_DIRS = [];
				COPY = True;
			elif line.strip() == "":
				COPY = False;
			elif COPY:
				TEST_DIRS.append(line.strip());
	return TEST_DIRS;

##############################
# FEATURE GENERATION SEGMENT
##############################
# mag_spec: generate magnitude spectrogram given audio
#	put features along rows
def gen_magspec(raw_audio, sr):

	#w_len = 0.02*sr; #(442,200)
	w_len = 0.006*sr; #(133, 669)
	#w_len = 0.0047*sr; #(104, 849)
	features =  librosa.stft(raw_audio,n_fft=np.int(w_len),hop_length=np.int(round(w_len/2)),window='hamming');
	#features = librosa.core.stft(raw_audio); #(173,1025)
	return features;

def gen_melspec(raw_audio,sr):
	n_fft = 1024
	hop_length = 512
	n_mels = 128

	# (133,669)
	#n_fft = np.int(0.006*sr);
	#hop_length = np.int(round(n_fft/2));
	#n_mels = 133;

	features = librosa.feature.melspectrogram(y=raw_audio,sr=sr,n_fft=n_fft,hop_length=hop_length,n_mels=n_mels);
	return features;

# generate modulation spec
def gen_modspec(raw_audio, sr, mtl_eng):
	modspec_out = mtl_eng.mymodspec_fn(raw_audio.tolist(),sr);
	features = np.asarray(modspec_out);

	return features;

# save_feat: save provided feature at provided output location
def save_feat(feat_in, out_path):
	feat = np.asarray(feat_in);
	feat = np.absolute(feat);
	with open(out_path,'w') as f: 
		np.savetxt(f,feat);

##############################
# MAIN CODE
##############################
### USER INPUT
# Define directory list filename
name = 'audio_all';
directory_list_filepath = name + '.dir_list.txt'; # which dataset to use
global_out = '/home/yichingl/mf/feat_out/'

# Define audio params
sr = 44100; # sampling rate
clip_duration = 2; # audio clips will split into segments of this duration (in seconds)

# Which features to use (only pick one bc they store to same directory)
processing_type = 'mod';
### EXTRACT FILENAMES CODE
# Construct pathways
train_paths = read_train_file(directory_list_filepath); # data input path
test_paths = read_test_file(directory_list_filepath); # data input path
print('Reading data from: '+ str(test_paths));

train_featout_paths = [line.replace("data","feat_out") for line in train_paths];
test_featout_paths = [line.replace("data","feat_out") for line in test_paths];
print('Exporting data to: '+ str(train_featout_paths));
print('Exporting data to: '+ str(test_featout_paths));

# Read filepaths_dict dict with keys 'train_filenames' 'train_labels' 'valid_filenames' 'valid_labels' 'test_filenames' 'test_labels'
#	and entries as list of filepaths (train/test_filenames)
#		and mats of labels (train/test_labels)
# TODO: see if this setup generates correctly, then add valid set before moving past file extraction

### START MATLAB ENGINE CODE
# This has to be in top-level file so it only has to start up once.. It takes a while
if processing_type == 'mod':
	print('Starting MATLAB engine...');
	mtl_eng = matlab.engine.start_matlab();
	mtl_eng.addpath(mtl_eng.genpath("matlab_toolboxes/"));
	print('Done importing MATLAB engine...');

# Read files from each of the training directories
train_ctl = '../ctl_files/' + name + '.' + processing_type + '.train.ctl';
with open(train_ctl, "w") as ctlout:
	n_dirs = len(train_paths);
	for n_dir in range(n_dirs):

		print('Reading data from: '+ str(train_paths[n_dir]));
		print('Writing data to: ' + str(train_featout_paths[n_dir]));

		filepaths_dict = read_filenames_set(train_paths[n_dir],train_featout_paths[n_dir]);
		in_filenames = filepaths_dict["in_filenames"];
		out_filenames = filepaths_dict["out_filenames"];

		n_files = len(in_filenames);
		for n_file in range(n_files):

			# select current file
			in_filename = in_filenames[n_file];
			out_filename = out_filenames[n_file];

			# check audio file's length in seconds (float) directly from file
			file_duration = librosa.get_duration(filename=in_filename);

			# load raw data using filename
			raw_data = load_audio(in_filename, sr);

			# calculate number of clips audio will be segmented into
			num_clips = int(np.floor(file_duration/clip_duration));
			print('num clips: '+ str(num_clips));

			# generate feature for for each input audio clip
			if num_clips != 0:
				start_ind = 0;
				clip_duration_dig = sr*clip_duration; # convert time in seconds to digital form
				for n in range(num_clips):
					end_ind = start_ind + clip_duration_dig - 1;

					# extract clip from raw data
					cur_clip = raw_data[start_ind:end_ind];

					if processing_type == 'spec':
						# generate magnitude spectrogram: has features along rows
						feat = gen_magspec(cur_clip, sr);

						# save generated feature to path
						out_clip_path = out_filename+"_"+str(n)+ ".spec";

					if processing_type == 'mel':
						# generate mel spectrogram: has features along rows
						feat = gen_melspec(cur_clip, sr);

						# save generated feature to path
						out_clip_path = out_filename+"_"+str(n)+ ".mel";

					if processing_type == 'mod':
						# generate mod spectrogram: has features along rows
						feat = gen_modspec(cur_clip, sr,mtl_eng);

						# save generated feature to path
						out_clip_path = out_filename+"_"+str(n)+ ".mod";

					print("	"+ out_clip_path);
					print("	"+ str(np.shape(feat)));
					ctlout.write(out_clip_path+"\n");
					
					save_feat(feat,out_clip_path);

					# update index
					start_ind = start_ind + clip_duration_dig;

			else:
				print("	provided clip is too short.. skipping...")
				
		print("Feature dimension: " +str(np.shape(feat)));


# Read files from each of the testing directories
test_ctl = '../ctl_files/' + name + '.'+processing_type+'.test.ctl';
with open(test_ctl, "w") as ctlout:
	n_dirs = len(test_paths);
	for n_dir in range(n_dirs):

		print('Reading data from: '+ str(test_paths[n_dir]));
		print('Writing data to: '+ str(test_featout_paths[n_dir]));

		filepaths_dict = read_filenames_set(test_paths[n_dir],test_featout_paths[n_dir]);
		in_filenames = filepaths_dict["in_filenames"];
		out_filenames = filepaths_dict["out_filenames"];

		n_files = len(in_filenames);
		for n_file in range(n_files):

			# select current file
			in_filename = in_filenames[n_file];
			print(in_filename);
			out_filename = out_filenames[n_file];

			# check audio file's length in seconds (float) directly from file
			file_duration = librosa.get_duration(filename=in_filename);

			# load raw data using filename
			raw_data = load_audio(in_filename, sr);

			# calculate number of clips audio will be segmented into
			num_clips = int(np.floor(file_duration/clip_duration));
			print('num clips: '+ str(num_clips));

			# generate feature for for each input audio clip
			if num_clips != 0:
				start_ind = 0;
				clip_duration_dig = sr*clip_duration; # convert time in seconds to digital form
				for n in range(num_clips):
					end_ind = start_ind + clip_duration_dig - 1;

					# extract clip from raw data
					cur_clip = raw_data[start_ind:end_ind];

					if processing_type == 'spec':
						# generate magnitude spectrogram: has features along rows
						feat = gen_magspec(cur_clip, sr);

						# save generated feature to path
						out_clip_path = out_filename+"_"+str(n)+ ".spec";

					if processing_type == 'mel':
						# generate mel spectrogram: has features along rows
						feat = gen_melspec(cur_clip, sr);

						# save generated feature to path
						out_clip_path = out_filename+"_"+str(n)+ ".mel";

					if processing_type == 'mod':
						# generate mod spectrogram: has features along rows
						feat = gen_modspec(cur_clip, sr,mtl_eng);

						# save generated feature to path
						out_clip_path = out_filename+"_"+str(n)+ ".mod";

					print("	"+ out_clip_path);
					print("	"+ str(np.shape(feat)));
					ctlout.write(out_clip_path+"\n");
					
					save_feat(feat,out_clip_path);

					# update index
					start_ind = start_ind + clip_duration_dig;

			else:
				print("	provided clip is too short.. skipping...")
				
		print("Feature dimension: " +str(np.shape(feat)));

