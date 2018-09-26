# for generating and saving graphics without needing display
# 	NOTE: this MUST be at the beginning of the file
# 	see: https://stackoverflow.com/questions/15713279/calling-pylab-savefig-without-display-in-ipython
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

from os import listdir # to get file names from a particular dir
import numpy as np
import librosa

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

# make_out_filenames: read all file names from given folder and return list of output filenames
def make_out_filenames(in_path,out_path):
	#filenames = [path + f for f in listdir(path) if (f[-4:] == ".wav")];  # include path
	filenames = [out_path + f[0:-4] for f in listdir(in_path) if (f[-4:] == ".wav")]; # don't include path
	return filenames;

# read_filenames_set: read in all files under train,test types in a dict with 
#	keys 'train' 'test'
#	and entries as lists of filenames
def read_filenames_set(datain_path, dataout_path,ttype,ctype):
	dataset = {};
	#bh = ['boat', 'heli'];
	#tvt = ['train',test];

	# read set (TODO: print to make sure feature appended correctly)
	filenames = read_filenames(datain_path+ctype+'/'+ttype+'/');
	out_filenames = make_out_filenames(datain_path+ctype+'/'+ttype+'/',dataout_path+ctype+'/'+ttype+'/');
		
	# store to output dict
	dataset[ttype+'_filenames'] = filenames;
	dataset[ttype+'_out_filenames'] = out_filenames;
	
	print('Reading set: (' +ttype+')');
	print('		'+str(len(filenames)) + ' files read');

	return dataset;

# load_audio: imports single audio file and return single raw data
def load_audio(filepath, sr, start_sec=0.0, duration_sec=None):
	# Import raw data
	try:
		raw_audio, sr = librosa.load(filepath, sr=sr, offset=start_sec, duration=duration_sec, mono=True);
	except:
		raise IOError('Give me an audio file which I can read!!');
	return raw_audio;


# mag_spec: generate magnitude spectrogram given audio
#	put features along rows
def gen_magspec(raw_audio, sr):

	#w_len = 0.02*sr; #442,200s
	#w_len = 0.006*sr;
	w_len = 0.0047*sr; #(104, 849)
	features =  librosa.stft(raw_audio,n_fft=np.int(w_len),hop_length=np.int(round(w_len/2)),window='hamming');
	#features = librosa.core.stft(raw_audio);
	return features;

def show_feat(feat, sr):
	print(np.shape(feat));
	fig = plt.figure();
	librosa.display.specshow(librosa.amplitude_to_db(feat, ref=np.max),x_axis='time',y_axis='linear');
	plt.title('Magnitude Spectrogram of Boat');
	plt.colorbar(format='%+2.0f dB');
	plt.tight_layout();
	plt.show();
	plt.close(fig);

# save_as_bin: save provided feature as a binary file at provided output
def save_as_bin(feat_in, out_path):
	feat = np.asarray(feat_in);
	feat = np.absolute(feat);
	num_valsi = np.size(feat);

	num_vals = np.asarray(num_valsi);

	with open(out_path, 'w') as f:
		num_vals.astype('int32').tofile(f,sep="",format="%i32");

	write_vals = np.transpose(feat).reshape(-1);

	with open(out_path, 'a') as f:
		write_vals.astype('float32').tofile(f,sep="",format="%f32");




##############################
# MAIN CODE
##############################

### USER INPUT
# Define file use params
dataset_name = 'audio_all'; # which dataset to use
navi_out = '../../../'; # path to get out of github directory
vis_out = 'graphics/'; # path to output visuals
global_out = '/home/yichingl/mf/feat_out/'

# Define audio params
sr = 44100; # sampling rate
clip_duration = 2; # audio clips will split into segments of this duration (in seconds)

# Define training params
max_epochs = 10; # maximum number of epochs for sgd
min_epochs = 2; # minimum number of epochs for sgd

### EXTRACT FILENAMES CODE
# Construct pathways
datain_path = navi_out + 'data/' + dataset_name + '/'; # data input path
featout_path = global_out + dataset_name + '/'; # data exported here

print('Reading data from: '+ datain_path);
print('Exporting data to: '+ featout_path);

# Read filepaths_dict dict with keys 'train_filenames' 'train_labels' 'valid_filenames' 'valid_labels' 'test_filenames' 'test_labels'
#	and entries as list of filepaths (train/test_filenames)
#		and mats of labels (train/test_labels)
# TODO: see if this setup generates correctly, then add valid set before moving past file extraction


# ### TRAINING CLASSIFIER CODE
# Use training set to build classifier
for ttype in ['train','test']:
	for ctype in ['boat','heli']:
		filepaths_dict = read_filenames_set(datain_path,featout_path,ttype,ctype);

		if(ttype=='train'):
			ttype_out = "../etc/"+ttype+"class."+ctype+".ctl";
		elif(ttype=='test'):
			ttype_out = ttype+"class."+ctype+".ctl";

		with open(ttype_out, "w") as ctlout:

			filenames = filepaths_dict[ttype+'_filenames'];
			out_filenames = filepaths_dict[ttype+'_out_filenames'];

			num_files = len(filenames);

			# Generate and save features
			for num_file in range(num_files):

				raw_audio_path = filenames[num_file];
				print("Reading from: "+raw_audio_path);
				out_path = out_filenames[num_file];
				

				# Check audio file's length in seconds (float) directly from file
				file_duration = librosa.get_duration(filename=raw_audio_path);
				#print('	file duration: '+ str(file_duration));

				# load raw data using filename
				raw_data = load_audio(raw_audio_path, sr);
				#print("File duration: "+str(np.shape(raw_data)));

				# Calculate number of clips audio will be segmented into
				num_clips = int(np.floor(file_duration/clip_duration));
				print('num clips: '+ str(num_clips));

				# generate feature for for each input audio clip]
				if num_clips != 0:
					start_ind = 0;
					clip_duration_dig = sr*clip_duration; # convert time in seconds to digital form
					for n in range(num_clips):
						end_ind = start_ind + clip_duration_dig - 1;

						# extract clip from raw data
						cur_clip = raw_data[start_ind:end_ind];
						#print("	clip duration: "+str(np.shape(cur_clip)));

						# generate magnitude spectrogram: has features along rows
						feat = gen_magspec(cur_clip, sr);

						# save generated feature to path
						out_clip_path = out_path+"_"+str(n)+ ".spec";
						print("	"+ out_clip_path);
						if(ttype=='train'):
							ctlout.write(out_clip_path+" "+str(1)+" "+str(np.shape(feat)[1])+"\n");
						elif(ttype=='test'):
							ctlout.write(out_clip_path+"\n");
						
						save_as_bin(feat,out_clip_path);

						# update index
						start_ind = start_ind + clip_duration_dig ;
				else:
					print("	provided clip is too short.. skipping...")
	print("Feature dimension: " +str(np.shape(feat)));

