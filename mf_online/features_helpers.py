import numpy as np
import librosa # to do signal processing
import audio_helpers as ah # to import audio
import anurag_cnn.feat_extractor as cnn # to import cnn
from anurag_cnn.feat_extractor import usegpu

# TODO: proofread this
# gen_features_multiple: generates a feature matrix using provided (multiple) audio filenames
def gen_features_multiple(filenames, labels, sr, clip_duration, mtl_eng):
	n_examples = len(filenames);
	
	# get dimensions (TODO: change this to something better?)
	n = 0;
	cur_data_dict = gen_features_single(filenames[n], labels[n], sr, clip_duration, mtl_eng);
	n_features = np.shape(cur_data_dict["features"])[1]; # number of cols

	# initialize features and labels arrays
	features_out = np.zeros((0,n_features));
	labels_out = np.zeros(0);

	for n in range(n_examples):

		# select current example
		cur_filename = filenames[n];
		cur_label = labels[n];

		# generate feature for current example
		cur_data_dict = gen_features_single(cur_filename, cur_label, sr, clip_duration, mtl_eng);
		if np.size(cur_data_dict["features"]) != 1:
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
def gen_features_single(filename, label, sr, clip_duration, mtl_eng):

	#print('Reading file: '+filename);

	# Check audio file's length in seconds (float) directly from file
	file_duration = librosa.get_duration(filename=filename);
	#print('	file duration: '+ str(file_duration));

	# load raw data using filename
	raw_data = ah.load_audio(filename, sr);

	# Calculate number of clips audio will be segmented into
	num_clips = int(np.floor(file_duration/clip_duration));
	print('	num clips: '+ str(num_clips));

	# generate feature for for each input audio clip]
	if num_clips != 0:
		start_ind = 0;
		clip_duration_dig = sr*clip_duration; # convert time in seconds to digital form
		for n in range(num_clips):
			end_ind = start_ind + clip_duration_dig - 1;

			# extract clip from raw data
			cur_clip = raw_data[start_ind:end_ind];

			# generate magnitude spectrogram: has features along rows
			magspec_feat = gen_magspec(cur_clip, sr);

			# generate cnn np array shape=(1,1024) and put along rows 
			cnn_feat = gen_cnn(cur_clip);
			cnn_feat = np.matlib.repmat(cnn_feat, np.shape(magspec_feat)[0], 1);

			# generate constant-Q, which has same time scale as magnitude spec
			constQ_feat = gen_constantQ(cur_clip, sr);

			# generate correlogram, which has same time scale as magnitude spec
			# TODO performance def decreased with addition of correlograms for miniset
			corr_feat = gen_correlogram(cur_clip,sr);

			# generate modulation spectrogram
			modspec_feat = gen_modspec(cur_clip,sr, mtl_eng);

			# combine generated feature types into one features matrix
			features_out_n = magspec_feat; # magnitude spectrogram
			features_out_n = np.hstack((features_out_n,cnn_feat)); # cnn
			features_out_n = np.hstack((features_out_n,constQ_feat)); # constant-Q
			features_out_n = np.hstack((features_out_n,corr_feat)); # correlogram
			# for mod spec
			features_out_n = np.repeat(features_out_n,2,0); #duplicate every row
			features_out_n = np.vstack((features_out_n[0,:],features_out_n));
			features_out_n = np.vstack((features_out_n,features_out_n[-1,:]));
			features_out_n = np.hstack((features_out_n,modspec_feat)); # modulation spectrogram

			# propagate label to all samples
			n_features = np.shape(features_out_n)[0];
			labels_out_n = np.multiply(label,np.ones(n_features, dtype=int));

			if n!=0:
				# if not first instance, append on items
				features_out = np.vstack((features_out,features_out_n));
				labels_out = np.concatenate((labels_out, labels_out_n));
			else:
				# if first instance, initialize
				features_out = features_out_n;
				labels_out = labels_out_n;

			# update index
			start_ind = start_ind + clip_duration_dig ;
	else:
		features_out = -11;
		labels_out = -11;

	# store to dict
	data_dict = {};
	data_dict["features"] = features_out;
	data_dict["labels"] = labels_out;

	return data_dict;

# generate correlogram, also generates in the same time scale as default
def gen_correlogram(raw_audio, sr):

    # Parameters
    acf_window = 10
    count_bins = 128 # Default for librosa, may want to change

    # Generate mel spectrograms
    spec = librosa.feature.melspectrogram(raw_audio, sr=sr, n_mels=count_bins)

    # Pad with zeros
    padding = np.zeros((count_bins, acf_window))
    spec_pad = np.concatenate((spec, padding), axis=1)

    # Generate chromagram
    #spec = librosa.feature.chroma_stft(data[0], sr=sampling_rate)
    
    # For each bin, calculate autocorrelation
    N = spec.shape[1]
    F = count_bins * min(N, acf_window)
    features = np.zeros((F, N))

    for i in range(N):
        obs = []
        for j in range(count_bins):
            #print("Len of autocorr: " + str(len(autocorr)))
            autocorr = librosa.core.autocorrelate(spec_pad[j, i:i+acf_window])[:N]
            obs.extend(autocorr)

        features[:, i] = np.transpose(obs)

    return np.transpose(features)

# generate modulation spec
def gen_modspec(raw_audio, sr, mtl_eng):
	features = np.asarray(mtl_eng.modspec_fn(raw_audio.tolist(),sr));

	return np.transpose(features);

# generate constant-Q
# 	put features along rows
# note: constantQ and stft generate in the same time scale, but might want to change to diff scale
def gen_constantQ(raw_audio, sr):
    features = np.abs(librosa.core.cqt(raw_audio, sr=sr));
    return np.transpose(features);
     
# gen_cnn: use a pre-trained CNN to extract features
def gen_cnn(raw_audio, use_gpu=True):
    
    features = cnn.main(raw_audio, use_gpu);
    
    return features;

# mag_spec: generate magnitude spectrogram given audio
#	put features along rows
def gen_magspec(raw_audio, sr):

	#w_len = 0.02*sr;
	#features =  librosa.stft(raw_audio,n_fft=w_len,hop_length=round(w_len/2),window='hamming');
	features = librosa.core.stft(raw_audio);
	return np.transpose(features);

def store_cnn(filenames_in, filenames_out, sr, clip_duration):
	n_clips = len(filenames_in);
	for n in range(n_clips):
		filename = filenames_in[n];
		featname = filenames_out[n];
		# Check audio file's length in seconds (float) directly from file
		file_duration = librosa.get_duration(filename=filename);
		#print('	file duration: '+ str(file_duration));

		# load raw data using filename
		raw_data = ah.load_audio(filename, sr);

		# Calculate number of clips audio will be segmented into
		num_clips = int(np.floor(file_duration/clip_duration));
		#print('	num clips: '+ str(num_clips));

		# generate feature for for each input audio clip
		start_ind = 0;
		clip_duration_dig = sr*clip_duration; # convert time in seconds to digital form
		for n in range(num_clips):
			end_ind = start_ind + clip_duration_dig - 1;

			# extract clip from raw data
			cur_clip = raw_data[start_ind:end_ind];

			# generate cnn np array shape=(1,1024) and put along rows 
			cnn_feat = gen_cnn(cur_clip);

			# combine generated feature types into one features matrix
			features_out_n = cnn_feat;


			if n!=0:
				# if not first instance, append on items
				features_out = np.vstack((features_out,features_out_n));
			else:
				# if first instance, initialize
				features_out = features_out_n;

			# update index
			start_ind = start_ind + clip_duration_dig ;

		# store output
		np.savetxt(featname, features_out, delimiter=',')   # X is an array
