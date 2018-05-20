import numpy as np
import librosa, librosa.display
import matplotlib.pyplot as plt


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

	#w_len = 0.02*sr;
	w_len = 0.005*sr;
	features =  librosa.stft(raw_audio,n_fft=np.int(w_len),hop_length=np.int(round(w_len/2)),window='hamming')
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
	print(np.shape(feat));
	print(np.mean(np.mean(feat)));
	print("0x{:016x}".format(88842)[2:]);

	num_vals = np.asarray(num_valsi);

	with open(out_path, 'w') as f:
		num_vals.astype('int32').tofile(f,sep="",format="%i32");

	#write_vals = np.transpose(feat).reshape(-1);
	#with open(out_path, 'w') as f:
	#write_vals.astype('float32').tofile(out_path,sep="",format="%s");
		#f.write(write_vals.tobytes())

	#print("number of vals " + str(num_vals) + ", " + str(np.binary_repr(num_vals, width=8)));

	#print("number of vals " + str(num_vals) + ", " + "0x{:016x}".format(num_vals)[2:]);
	#write_vals = np.append(np.array([num_vals]), np.transpose(feat).reshape(-1));
	write_vals = np.transpose(feat).reshape(-1);
	#write_bin = np.asarray(["0x{:016x}".format(x)[2:] for x in write_vals]);

	#write_bin = [[format(x,'08b') for x in y] for y in write_vals] for 2d lists#



	with open(out_path, 'a') as f:
		write_vals.astype('float32').tofile(f,sep="",format="%f32");

	# out-test
	#out_test_path = '/home/yichingl/mf/data/miniset/heli/test/feats/test.spec';

	#test_num_vals = 5;

	#with open(out_path, 'w') as f:
	#	test_num_vals.astype('int32').tofile(f,sep="",format="%i32");

	#with open(out_path, 'a') as f:
	#	write_vals.astype('float32').tofile(f,sep="",format="%f32");





	#write_vals.tofile(out_path,sep="");

raw_audio_path = '../../../data/miniset/heli/test/72qNna11Znw_051.wav';
out_path = '../../../data/miniset/heli/test/feats/72qNna11Znw_051.spec';
sr = 44100;

raw_audio = load_audio(raw_audio_path,sr, start_sec=0.0,duration_sec=2);
feat = gen_magspec(raw_audio,sr);
show_feat(feat,sr);
#save_as_bin(feat,out_path);