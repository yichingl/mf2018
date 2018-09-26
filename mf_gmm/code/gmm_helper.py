import numpy as np;
from sklearn import mixture;
import bob.learn.em

# reads ctl file for specifiec spec type
def read_ctlname(line,feat_type):
	line = line.rstrip('\n');
	ctl_name = line.replace('spec',feat_type);
	return ctl_name

# reads ctl file and returns numpy array containing data [n_samples, n_feats]
def read_specgram_features(ctl_path, n_feats, feat_type):
	feat_all = np.zeros((n_feats,0));
	with open(ctl_path) as f: 
		for line in f:
			ctl_name = read_ctlname(line, feat_type);
			print('Reading (train gmm): ' + ctl_name);

			# read data
			read_data = np.absolute(np.loadtxt(ctl_name)); 
		feat_all = np.hstack((feat_all,read_data));
	return np.transpose(feat_all) #[n_feats, n_samples] -> [n_samples, n_feats]

# read ctl file and return spec magnitude data
def read_specgram_feature(ctl_name):
	print('Reading (supervector-train): ' + ctl_name);
	read_data = np.absolute(np.loadtxt(ctl_name)); 
	return np.transpose(read_data);

# run k-means on data and return means to use for GMM MLE algorithm
def get_kmeans_means(feat_all, n_clusters, n_feats):

	# parameters
	max_iterations = 200;
	convergence_threshold = 1e-5;

	# initialize kmeans machine and trainer
	kmeans_machine = bob.learn.em.KMeansMachine(n_clusters, n_feats);
	kmeans_trainer = bob.learn.em.KMeansTrainer();

	# train the KMeansMachine
	bob.learn.em.train(kmeans_trainer, kmeans_machine, feat_all, max_iterations=max_iterations,convergence_threshold=convergence_threshold);

	# extract means
	kmeans_means = kmeans_machine.means;

	return kmeans_means;

# trains GMM and returns the machine
def train_gmm(feat_all, n_clusters, n_feats):

	# parameters
	max_iterations = 200;
	convergence_threshold = 1e-5;

	# Create a kmeans model (machine) with n_clusters and dimensionality equal to n_feats
	gmm_ubm = bob.learn.em.GMMMachine(n_clusters,n_feats); 
	# Use MLE trainer to train the GMM 
	# (True,True,True means update means/variances/weights at each iteration)
	gmm_trainer = bob.learn.em.ML_GMMTrainer(True, True, True);
	# Train using MLE trainer
	bob.learn.em.train(gmm_trainer, gmm_ubm, feat_all, max_iterations=max_iterations,convergence_threshold=convergence_threshold);

	return (gmm_ubm)

# perform MAP adaptation and return new means
def create_supervect(feat_read, n_clusters, n_feats, gmm_ubm):

	# parameters
	max_iterations = 200;
	convergence_threshold = 1e-5;

	# Create a new GMM model (machine) to adapt using MAP
	adapted_gmm = bob.learn.em.GMMMachine(n_clusters,n_feats);

	# Create the MAP trainer
	gmm_trainer = bob.learn.em.MAP_GMMTrainer(gmm_ubm, relevance_factor=4);

	# Train using MAP
	bob.learn.em.train(gmm_trainer, adapted_gmm, feat_read, max_iterations=max_iterations,convergence_threshold=convergence_threshold);

	supervector = adapted_gmm.mean_supervector; # (n_means, n_inputs) 

	return supervector

# # perform MAP adaptation and return new means
# def create_supervect(feat_read, n_clusters, n_feats, init_means):

# 	# parameters
# 	max_iterations = 200;
# 	convergence_threshold = 1e-5;

# 	# Create a prior GMM model (machine) with n_clusters and dimensionality equal to n_feats
# 	gmm_ubm = bob.learn.em.GMMMachine(n_clusters,n_feats);
# 	# Set the means
# 	gmm_ubm.means = init_means;

# 	# Create a new GMM model (machine) to adapt using MAP
# 	adapted_gmm = bob.learn.em.GMMMachine(n_clusters,n_feats);

# 	# Create the MAP trainer
# 	gmm_trainer = bob.learn.em.MAP_GMMTrainer(gmm_ubm, relevance_factor=4);

# 	# Train using MAP
# 	bob.learn.em.train(gmm_trainer, adapted_gmm, feat_read, max_iterations=max_iterations,convergence_threshold=convergence_threshold);

# 	supervector = adapted_gmm.means; # (n_means, n_inputs) 
# 	supervector = supervector.flatten('F'); # flatten column major

# 	return supervector
