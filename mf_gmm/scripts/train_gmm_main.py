
import numpy as np;
from sklearn import mixture;

# read binary file located at path
def read_bin(path, nvects):
	num_vals = int(np.fromfile(path,dtype=np.int32,count=1));
	feat_vect = np.fromfile(path,dtype=np.float32,count=-1);
	feat_vect = feat_vect[1:];

	# check if input and num_vals match?
	if(num_vals != np.size(feat_vect)):
		print("values not equal");

	feat_out = np.reshape(feat_vect,(num_vals/nvects,nvects));

	return np.transpose(feat_out);




##########################################################
# Main: read list of files from .ctl, train gmm and store
##########################################################
# don't use start and end for now
dataset_name = 'audio_all'
gmm_weights_path = dataset_name+'.gmm.weights.csv';
gmm_means_path = dataset_name+'.gmm.means.csv';

# initialize gmm
n_clusters = 5;
gmm = mixture.GaussianMixture(n_components=n_clusters);

ctlset_names = ['trainclass.boat.ctl','trainclass.heli.ctl'];

ctlset_path = '../etc';

feat_all = np.zeros((0,104));

# loop through control files
for ctlset_name in ctlset_names:
	with open(ctlset_path+'/'+ctlset_name) as f: 

		# extract file name, and start and end vectors
		for line in f:
			line = line.rstrip('\n');
			txt = line.split(" ");
			ctl_name = txt[0];
			start_n = int(txt[1]);
			end_n = int(txt[2]);
			print('Reading: ' + ctl_name);

			# read data from control file
			feat_all = np.vstack((feat_all,read_bin(ctl_name, end_n)));
gmm.fit(feat_all);
print(np.asarray(gmm.means_));
# save gmm model
with open(gmm_weights_path, 'w') as f:
	np.savetxt(f,np.asarray(gmm.weights_));
with open(gmm_means_path, 'w') as f:
	np.savetxt(f,np.asarray(gmm.means_));
