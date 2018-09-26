
import numpy as np;
from sklearn import mixture;

# read binary file located at path
def read_bin1(path, nfeats):
	num_vals = int(np.fromfile(path,dtype=np.int32,count=1));
	feat_vect = np.fromfile(path,dtype=np.float32,count=-1);
	feat_vect = feat_vect[1:];

	# check if input and num_vals match?
	if(num_vals != np.size(feat_vect)):
		print("values not equal");

	feat_out = np.reshape(feat_vect,(nfeats,num_vals/nfeats));

	return np.transpose(feat_out);

# read binary file located at path
def read_bin2(path, nvects):
	num_vals = int(np.fromfile(path,dtype=np.int32,count=1));
	feat_vect = np.fromfile(path,dtype=np.float32,count=-1);
	feat_vect = feat_vect[1:];

	# check if input and num_vals match?
	if(num_vals != np.size(feat_vect)):
		print("values not equal");

	feat_out = np.reshape(feat_vect,(num_vals/nvects,nvects));

	return np.transpose(feat_out);


##########################################################
# Main: read list of files from .ctl, create supervects and store
##########################################################
# don't use start and end for now

### USER INPUT
# Define file use params
n_clusters = 5;

dataset_name = 'audio_all'; # which dataset to use
navi_out = '../../../'; # path to get out of github directory
global_out = '/home/yichingl/mf/supervect_out/'

gmm_weights_path = dataset_name+'.gmm.weights.csv';
gmm_means_path = dataset_name+'.gmm.means.csv';

# Construct pathways
featout_path = global_out + dataset_name + '/'; # data exported here

# load gmm
init_weights = np.loadtxt(gmm_weights_path);
init_means = np.loadtxt(gmm_means_path);

ctlset_names_test = ['testclass.boat.ctl','testclass.heli.ctl'];
ctlset_names_train = ['trainclass.boat.ctl','trainclass.heli.ctl'];

ctlset_path = '../etc';

nfeats = 104;


# loop through control files
for ctlset_name in ctlset_names_test:
	with open(ctlset_path+'/'+ctlset_name) as f: 

		# extract file name, and start and end vectors
		for line in f:
			ctl_name = line.rstrip('\n');
			print('Reading: ' + ctl_name);

			# read data from control file into 2D matrix form
			feat_read = read_bin1(ctl_name,nfeats);

			# update supervector gmm once
			gmm = mixture.GaussianMixture(n_components=n_clusters, max_iter=1, weights_init=init_weights, means_init=init_means);
			gmm.fit(feat_read);

			# save means
			supervect_out_path = ctl_name.replace('feat_out','supervect_out');
			with open(supervect_out_path, 'w') as f:
				np.savetxt(f,np.asarray(gmm.means_));


# loop through control files
for ctlset_name in ctlset_names_train:
	with open(ctlset_path+'/'+ctlset_name) as f: 

		# extract file name, and start and end vectors
		for line in f:
			line = line.rstrip('\n');
			txt = line.split(" ");
			ctl_name = txt[0];
			start_n = int(txt[1]);
			end_n = int(txt[2]);
			print('Reading: ' + ctl_name);

			# read data from control file into 2D matrix form
			feat_read = read_bin2(ctl_name,end_n);

			# update supervector gmm once
			gmm = mixture.GaussianMixture(n_components=n_clusters, max_iter=1, weights_init=init_weights, means_init=init_means);
			gmm.fit(feat_read);

			# save means
			supervect_out_path = ctl_name.replace('feat_out','supervect_out');
			with open(supervect_out_path, 'w') as f:
				np.savetxt(f,np.asarray(gmm.means_));