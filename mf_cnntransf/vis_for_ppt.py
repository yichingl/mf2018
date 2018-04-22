import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from os import listdir # to get file names from a particular dir
import numpy as np
import audio_helpers as ah
import features_helpers as fh

import librosa, librosa.display # to plot spectrograms

##############################
# MAIN CODE
##############################

### USER INPUT
# Define file use params
dataset_name = 'short'; # which dataset to use
navi_out = '../../'; # path to get out of github directory
vis_out = 'graphics/'; # path to output visuals

# Define audio params
sr = 44100; # sampling rate
clip_duration = 2; # audio clips will split into segments of this duration (in seconds)

# Define training params
max_epochs = 2; # maximum number of epochs for sgd
min_epochs = 1; # minimum number of epochs for sgd


### EXTRACT FILENAMES CODE
# Construct pathways
datain_path = navi_out + 'data/' + dataset_name + '/'; # data input path
dataout_path = navi_out + 'data_out/' + dataset_name + '/'; # data exported here
vis_path = vis_out + dataset_name+'/'; # graphics exported here within repo

print('Reading data from: '+ datain_path);
print('Exporting data to: '+ dataout_path);

filepath_boat = datain_path+'boat/boat9_short.wav';
raw_data_boat = ah.load_audio(filepath_boat,sr=sr,start_sec=0,duration_sec=clip_duration);
filepath_heli = datain_path+'heli/heli6_short.wav';
raw_data_heli = ah.load_audio(filepath_heli,sr=sr,start_sec=0,duration_sec=clip_duration);

print(len(raw_data_boat));
# Generate constantQ spectrogram image
# feat = np.transpose(fh.gen_constantQ(raw_data_boat,sr));
# fig = plt.figure();
# librosa.display.specshow(librosa.amplitude_to_db(feat, ref=np.max),x_axis='time');
# plt.title('Boat ConstantQ Spectrogram');
# plt.colorbar(format='%+2.0f dB');
# plt.tight_layout();
# plt.savefig(vis_path+'constantq__boat9.png');
# plt.close(fig);
# feat = np.transpose(fh.gen_constantQ(raw_data_heli,sr));
# fig = plt.figure();
# librosa.display.specshow(librosa.amplitude_to_db(feat, ref=np.max),x_axis='time');
# plt.title('Helicopter ConstantQ Spectrogram');
# plt.colorbar(format='%+2.0f dB');
# plt.tight_layout();
# plt.savefig(vis_path+'constantq__heli6.png');
# plt.close(fig);
feat = librosa.amplitude_to_db(librosa.cqt(raw_data_boat,sr), ref=np.max);
fig = plt.figure();
librosa.display.specshow(feat, y_axis='cqt_hz',x_axis='time',sr=sr);
plt.title('Boat ConstantQ Spectrogram');
plt.colorbar(format='%+2.0f dB');
plt.tight_layout();
plt.savefig(vis_path+'constantq_boat9.png');
plt.close(fig);
feat = librosa.amplitude_to_db(librosa.cqt(raw_data_heli,sr), ref=np.max);
fig = plt.figure();
librosa.display.specshow(feat,x_axis='time', y_axis='cqt_hz',sr=sr);
plt.title('Helicopter ConstantQ Spectrogram');
plt.colorbar(format='%+2.0f dB');
plt.tight_layout();
plt.savefig(vis_path+'constantq_heli6.png');
plt.close(fig);

# Generate correlogram image
feat = fh.gen_correlogram(raw_data_boat,sr);
fig = plt.figure();
librosa.display.specshow(librosa.amplitude_to_db(feat, ref=np.max),x_axis='time',y_axis='mel',sr=sr);
plt.title('Boat Correlogram');
plt.colorbar(format='%+2.0f dB');
plt.tight_layout();
plt.savefig(vis_path+'correlogram_boat9.png');
plt.close(fig);
feat = fh.gen_correlogram(raw_data_heli,sr);
fig = plt.figure();
librosa.display.specshow(librosa.amplitude_to_db(feat, ref=np.max),x_axis='time',y_axis='mel',sr=sr);
plt.title('Heli Correlogram');
plt.colorbar(format='%+2.0f dB');
plt.tight_layout();
plt.savefig(vis_path+'correlogram_heli6.png');
plt.close(fig);

# Generate correlogram image
feat = fh.gen_mel(raw_data_boat,sr);
fig = plt.figure();
librosa.display.specshow(librosa.amplitude_to_db(feat, ref=np.max),x_axis='time',y_axis='mel',sr=sr);
plt.title('Boat Melspectrogram');
plt.colorbar(format='%+2.0f dB');
plt.tight_layout();
plt.savefig(vis_path+'melspecgram_boat9.png');
plt.close(fig);
feat =fh.gen_mel(raw_data_heli,sr);
fig = plt.figure();
librosa.display.specshow(librosa.amplitude_to_db(feat, ref=np.max),x_axis='time',y_axis='mel',sr=sr);
plt.title('Heli Melspectrogram');
plt.colorbar(format='%+2.0f dB');
plt.tight_layout();
plt.savefig(vis_path+'melspecgram_heli6.png');
plt.close(fig);
