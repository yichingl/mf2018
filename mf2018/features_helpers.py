import numpy as np
import librosa # to do signal processing
import audio_helpers as ah # to import audio
import anurag_cnn.feat_extractor as cnn # to import cnn
from anurag_cnn.feat_extractor import usegpu

# TODO: proofread this
# gen_features_multiple: generates a feature matrix using provided (multiple) audio filenames
def gen_features_multiple(filenames, labels, sr, clip_duration):
	n_examples = len(filenames);
	
	# get dimensions (TODO: change this to something better?)
	n = 0;
	cur_data_dict = gen_features_single(filenames[n], labels[n], sr, clip_duration);
	n_features = np.shape(cur_data_dict["features"])[1]; # number of cols

	# initialize features and labels arrays
	features_out = np.zeros((0,n_features));
	labels_out = np.zeros(0);

	for n in range(n_examples):

		# select current example
		cur_filename = filenames[n];
		cur_label = labels[n];

		# generate feature for current example
		cur_data_dict = gen_features_single(cur_filename, cur_label, sr, clip_duration);
		features_out = np.vstack((features_out, cur_data_dict["features"]));
		labels_out = np.concatenate((labels_out, cur_data_dict["labels"]));

	data_dict = {};
	data_dict["features"] = features_out;
	data_dict["labels"] = labels_out;

	return data_dict;


# gen_features_single: generates a feature matrix using provided (single) audio filename
#	inputs:
#		filename
#		label: single label to be propagated for all output examples
#		sr: sampling rate (use to convert between analog and digital)
#		clip_duration: duration to segment data into (used for various feature extraction methods)
#	output:
#		data_dict: holds "features" and "labels"
#			"features" will have each feature vector along rows
def gen_features_single(filename, label, sr, clip_duration):

	# load raw data using filename
	raw_data = ah.load_audio(filename, sr);

	# generate magnitude spectrogram and put features along rows
	mag_spec = np.transpose(gen_mag_spec(raw_data, sr));

	# combine generated features into one features matrix
	features_out = mag_spec;
	
	# propagate label to all samples
	n_features = np.shape(features_out)[0];
	labels_out = np.multiply(label,np.ones(n_features, dtype=int));

	# store to dict
	data_dict = {};
	data_dict["features"] = features_out;
	data_dict["labels"] = labels_out;

	return data_dict;

# gen_cnn: use a pre-trained CNN to extract features
def gen_cnn(raw_audio, use_gpu=True):
    
    features = cnn.main(raw_audio[0], use_gpu)
    
    return features
# mag_spec: generate magnitude spectrogram given audio
def gen_mag_spec(raw_audio, sr):

	return librosa.stft(raw_audio);
	w_len = 0.02*sr;
	return librosa.stft(raw_audio,n_fft=w_len,hop_length=round(w_len/2),window='hamming');