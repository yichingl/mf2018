# for generating and saving graphics without needing display
# 	NOTE: this MUST be at the beginning of the file
# 	see: https://stackoverflow.com/questions/15713279/calling-pylab-savefig-without-display-in-ipython
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from os import listdir # to get file names from a particular dir
import numpy as np
import features_helpers as fh
import audio_helpers as ah # to import audio
import librosa, librosa.display # to plot spectrograms


# Define audio params
sr = 44100; # sampling rate
clip_duration = 2; # audio clips will split into segments of this duration (in seconds)

# Define file use params
dataset_name = 'miniset'; # which dataset to use
navi_out = '../../'; # path to get out of github directory
vis_out = 'graphics/'; # path to output visuals

# Construct pathways
datain_path = navi_out + 'data/' + dataset_name + '/'; # data input path
vis_path = vis_out + dataset_name+'/'; # graphics exported here within repo

filepath_boat = datain_path+'boat/test/v55hojInS9c_011.wav';
raw_data_boat = ah.load_audio(filepath_boat,sr=sr,start_sec=0,duration_sec=clip_duration);
filepath_heli = datain_path+'heli/test/72qNna11Znw_045.wav';
raw_data_heli = ah.load_audio(filepath_heli,sr=sr,start_sec=0,duration_sec=clip_duration);

# Generate magnitude spectrogram image
feat = np.transpose(fh.gen_magspec(raw_data_boat,sr));
fig = plt.figure();
librosa.display.specshow(librosa.amplitude_to_db(feat, ref=np.max),x_axis='time');
plt.title('Magnitude Spectrogram of Boat');
plt.colorbar(format='%+2.0f dB');
plt.tight_layout();
plt.savefig(vis_path+'magspec__boat_v55hojInS9c_011.png');
plt.close(fig);
feat = np.transpose(fh.gen_magspec(raw_data_heli,sr));
fig = plt.figure();
librosa.display.specshow(librosa.amplitude_to_db(feat, ref=np.max),x_axis='time');
plt.title('Magnitude Spectrogram of Heli');
plt.colorbar(format='%+2.0f dB');
plt.tight_layout();
plt.savefig(vis_path+'magspec__heli_72qNna11Znw_045.png');
plt.close(fig);

# Generate constantQ spectrogram image
feat = np.transpose(fh.gen_constantQ(raw_data_boat,sr));
fig = plt.figure();
librosa.display.specshow(librosa.amplitude_to_db(feat, ref=np.max),x_axis='time');
plt.title('ConstantQ Spectrogram of Boat');
plt.colorbar(format='%+2.0f dB');
plt.tight_layout();
plt.savefig(vis_path+'constantq__boat_v55hojInS9c_011.png');
plt.close(fig);
feat = np.transpose(fh.gen_constantQ(raw_data_heli,sr));
fig = plt.figure();
librosa.display.specshow(librosa.amplitude_to_db(feat, ref=np.max),x_axis='time');
plt.title('ConstantQ Spectrogram of Heli');
plt.colorbar(format='%+2.0f dB');
plt.tight_layout();
plt.savefig(vis_path+'constantq__heli_72qNna11Znw_045.png');
plt.close(fig);

# Generate correlogram image
feat = np.transpose(fh.gen_correlogram(raw_data_boat,sr));
fig = plt.figure();
librosa.display.specshow(librosa.amplitude_to_db(feat, ref=np.max),x_axis='time');
plt.title('Correlogram of Boat');
plt.colorbar(format='%+2.0f dB');
plt.tight_layout();
plt.savefig(vis_path+'correlogram__boat_v55hojInS9c_011.png');
plt.close(fig);
feat = np.transpose(fh.gen_correlogram(raw_data_heli,sr));
fig = plt.figure();
librosa.display.specshow(librosa.amplitude_to_db(feat, ref=np.max),x_axis='time');
plt.title('Correlogram of Heli');
plt.colorbar(format='%+2.0f dB');
plt.tight_layout();
plt.savefig(vis_path+'correlogram__heli_72qNna11Znw_045.png');
plt.close(fig);

#Generate CNN
feat = fh.gen_cnn(raw_data_boat);
metrics_path = vis_path+"cnn__boat_v55hojInS9c_011.csv";
np.savetxt(metrics_path, feat, delimiter=",",fmt='%f');
feat = fh.gen_cnn(raw_data_heli);
metrics_path = vis_path+"cnn__heli_72qNna11Znw_045.csv";
np.savetxt(metrics_path, feat, delimiter=",",fmt='%f');