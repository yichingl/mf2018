
from os import listdir # to get file names from a particular dir

# input ctl name to save
folders = ['/home/yichingl/mf/data/minitest/boat','/home/yichingl/mf/data/minitest/heli'];

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
# MAIN CODE
##############################
### USER INPUT
# Define directory list filename
global_out = '/home/yichingl/mf/feat_out/'
for directory in folders:

	# Define audio params
	sr = 44100; # sampling rate
	clip_duration = 2; # audio clips will split into segments of this duration (in seconds)

	# Which features to use (only pick one bc they store to same directory)
	processing_type = 'mel';
	### EXTRACT FILENAMES CODE
	# Construct pathways
	train_paths = read_train_file(directory_list_filepath); # data input path
	test_paths = read_test_file(directory_list_filepath); # data input path
	print('Reading data from: '+ str(test_paths));

	train_featout_paths = [line.replace("data","feat_out") for line in train_paths];
	test_featout_paths = [line.replace("data","feat_out") for line in test_paths];
	print('Exporting data to: '+ str(train_featout_paths));

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
