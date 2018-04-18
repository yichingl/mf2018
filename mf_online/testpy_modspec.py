import numpy as np
import librosa # to do signal processing
import matlab.engine

print('Starting MATLAB engine...');
mtl_eng = matlab.engine.start_matlab();
mtl_eng.addpath(mtl_eng.genpath("matlab_toolboxes/"));
print('Done importing MATLAB engine...');

#import audio
print('Starting code...');
sr = 44100;
raw_audio,sr = librosa.load("../../data/miniset/boat/train/ECyqCpZoXWg_009.wav",sr=sr);
raw_audio = raw_audio[0:2*sr-1];
modspecgram = np.asarray(mtl_eng.modspec_fn(raw_audio.tolist(),sr));
print('Completed code..');
print(np.shape(modspecgram));