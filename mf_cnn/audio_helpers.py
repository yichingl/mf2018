import numpy as np
import librosa # to load audio


# load_audio: imports single audio file and return single raw data
def load_audio(filepath, sr, start_sec=0.0, duration_sec=None):
	# Import raw data
	try:
		raw_audio, sr = librosa.load(filepath, sr=sr, offset=start_sec, duration=duration_sec, mono=True);
	except:
		raise IOError('Give me an audio file which I can read!!');
	return raw_audio;