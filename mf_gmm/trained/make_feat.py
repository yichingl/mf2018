import numpy as np;
import librosa; # to load audio and do signal processing


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

	w_len = 0.02*sr;
	features =  librosa.stft(raw_audio,n_fft=np.int(w_len),hop_length=np.int(round(w_len/2)),window='hamming')
	#features = librosa.core.stft(raw_audio);
	return features;

# save_as_bin: save provided feature as a binary file at provided output
def save_as_bin(feat_in, out_path):
	feat = np.asarray(feat_in);
	num_vals = np.size(feat);
	write_vals = np.append(num_vals, np.transpose(feat).reshape(-1));
	write_vals.tofile(out_path,sep="",format="%int8");

raw_audio_path = '../../data/mini/heli/test/heli8_short.wav';
out_path = '../../data/mini/heli/test/feats/heli8_specgram_bin.spec';
sr = 44100;

raw_audio = load_audio(raw_audio_path,sr, start_sec=0.0,duration_sec=2);
feat = gen_magspec(raw_audio,sr);
print('hi');
print(np.shape(feat));
save_as_bin(feat,out_path);