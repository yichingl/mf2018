import numpy as np;
from sklearn import mixture;
from sklearn import svm;
from sklearn.linear_model import SGDClassifier;
import random; # for shuffling

TEST = 0;

BOAT_N = 0;
HELI_N = 1;
##############################
# TRAIN GMM HELPER CODE
##############################

##############################
# MAIN CODE
##############################

# USER INPUT
# For gmm
group_name='sel';
ctl_group_names = ['sel.etc','sel.boatA','sel.heliA'];
ctl_train_names = ['sel.boatA','sel.heliA']
ctl_test_names = ['sel.boatA','sel.heliA']
# For svm
ctl_train_boat = 'sel.boatA.train.ctl';
ctl_train_heli = 'sel.heliA.train.ctl';
ctl_test_boat = 'sel.boatA.test.ctl';
ctl_test_heli = 'sel.heliA.test.ctl';

if TEST:
	print('TESTING');
	group_name = 'test';
	ctl_group_names = ['test.boat','test.heli'];
	ctl_train_names = ['test.boat','test.heli']
	ctl_train_boat = 'test.boat.train.ctl';
	ctl_train_heli = 'test.heli.train.ctl';
	ctl_test_boat = 'test.boat.test.ctl';
	ctl_test_heli = 'test.heli.test.ctl';

feat_out = '/home/yichingl/mf/feat_out/';
supervect_out = '/home/yichingl/mf/supervect_out/';

feat_type = 'spec';

# Define tasks
TRAIN_GMM = 1;
MAKE_SUPERVECTS = 1;
CLASSIFY = 1;

# Define audio params
sr = 44100; # sampling rate
clip_duration = 2; # audio clips will split into segments of this duration (in seconds)

# Define gmm params
n_clusters = 64;
max_iter = 5;
covariance_type = 'diag';

# Define svm params
svm_type = 'full';
C = .85; #default 1
kernel = 'linear';#rbf''poly' 'linear'
batch_size_opts = [500,1000,2500];
n_iters_opts = [100,500,1000];
if svm_type == 'full':
	batch_size_opts = ['full'];
	n_iters_opts = ['full'];





# hardcode (specgram)
if feat_type == 'spec':
	#n_feats = 442;
	#n_timeframes = 200;

	#n_feats = 104;
	#n_timeframes = 849;

	#n_feats = 133;
	#n_timeframes = 669;

	n_feats = 1025;
	n_timeframes = 173;

elif feat_type == 'mel':
	n_feats = 128;
	n_timeframes = 173;
elif feat_type == 'mod':
	n_feats = 128;
	n_timeframes = 320;

##############################
# TRAIN GMM CODE
##############################

# Read train ctl file
if TRAIN_GMM:
	# User input

	for ctl_group_name in ctl_group_names:
		# don't use start and end for now
		gmm_weights_path = '../trained_gmms/'+ctl_group_name+'.'+feat_type+'.cl'+str(n_clusters)+'.gmm.weights.csv';
		gmm_means_path = '../trained_gmms/'+ctl_group_name+'.'+feat_type+'.cl'+str(n_clusters)+'.gmm.means.csv';
		ctl_path = '../ctl_files/' + ctl_group_name +'.train.ctl';

		# initialize gmm
		gmm = mixture.GaussianMixture(n_components=n_clusters, max_iter=max_iter, covariance_type=covariance_type);

		# hardcode
		feat_all = np.zeros((0,n_feats));

		# loop through control file
		with open(ctl_path) as f: 

			# extract file name, and start and end vectors
			for ctl_name in f:
				ctl_name = ctl_name.rstrip('\n');
				ctl_name = ctl_name.replace('spec',feat_type);
				print('Reading (train gmm): ' + ctl_name);

				# read data
				x = np.transpose(np.loadtxt(ctl_name));

				# normalize data
				x = x/(np.linalg.norm(x, ord=1, axis=1, keepdims=True)+np.finfo(float).eps);
				print(np.shape(x));

				# read data from file and stack
				feat_all = np.vstack((feat_all,x));

		gmm.fit(feat_all);
		# save gmm model
		with open(gmm_weights_path, 'w') as f:
			np.savetxt(f,np.asarray(gmm.weights_));
		with open(gmm_means_path, 'w') as f:
			np.savetxt(f,np.asarray(gmm.means_));

##############################
# MAKE SUPERVECTOR CODE
##############################
if MAKE_SUPERVECTS:

	for ctl_group_name_supervect in ctl_train_names:
		# don't use start and end for now
		gmm_weights_path = '../trained_gmms/'+ctl_group_name_supervect+'.'+feat_type+'.cl'+str(n_clusters)+'.gmm.weights.csv';
		gmm_means_path = '../trained_gmms/'+ctl_group_name_supervect+'.'+feat_type+'.cl'+str(n_clusters)+'.gmm.means.csv';

		try:
			init_weights = np.stack((init_weights,np.loadtxt(gmm_weights_path)),axis=0);
			init_means = np.stack((init_means,np.loadtxt(gmm_means_path)),axis=0);
		except NameError:
			# load gmm
			init_weights = np.loadtxt(gmm_weights_path);
			init_means = np.loadtxt(gmm_means_path);

	# agglomarate and normalize
	init_weights = np.mean(init_weights,axis=0);
	init_weights = init_weights/(np.linalg.norm(init_weights, ord=1)+np.finfo(float).eps);

	init_means = np.mean(init_means,axis=0);
	init_means = init_means/(np.linalg.norm(init_means, ord=1, axis=1, keepdims=True)+np.finfo(float).eps);

	for ctl_group_name_supervect in ctl_test_names:
		ctl_path = '../ctl_files/' + ctl_group_name_supervect +'.train.ctl';
		# loop through control file
		with open(ctl_path) as f: 
			# extract file name, and start and end vectors
			for line in f:
				ctl_name = line.rstrip('\n');
				ctl_name = ctl_name.replace('spec',feat_type);
				print('Reading (supervector-train): ' + ctl_name);

				# read data from control file into 2D matrix form
				feat_read = np.transpose(np.loadtxt(ctl_name));

				# update supervector gmm once
				gmm = mixture.GaussianMixture(n_components=n_clusters, max_iter=1, covariance_type=covariance_type, weights_init=init_weights, means_init=init_means);
				gmm.fit(feat_read);

				# save means
				supervect_out_path = ctl_name.replace('feat_out','supervect_out');
				with open(supervect_out_path, 'w') as f:
					np.savetxt(f,np.asarray(gmm.means_));

	ctl_path = '../ctl_files/' + 'sel' +'.test.ctl';
	# loop through control file
	with open(ctl_path) as f: 
		# extract file name, and start and end vectors
		for line in f:
			ctl_name = line.rstrip('\n');
			ctl_name = ctl_name.replace('spec',feat_type);
			print('Reading (supervector-test): ' + ctl_name);

			# read data from control file into 2D matrix form
			feat_read = np.transpose(np.loadtxt(ctl_name));
			print(np.shape(feat_read));

			# update supervector gmm once
			gmm = mixture.GaussianMixture(n_components=n_clusters, max_iter=1, covariance_type=covariance_type, weights_init=init_weights, means_init=init_means);
			gmm.fit(feat_read);

			# save means
			supervect_out_path = ctl_name.replace('feat_out','supervect_out');
			with open(supervect_out_path, 'w') as f:
				np.savetxt(f,np.asarray(gmm.means_));


for batch_size in batch_size_opts:
	for n_iters in n_iters_opts:
		print('\n\n\nn_clusters: '+str(n_clusters)+', batch_size: '+str(batch_size)+', n_iters: '+str(n_iters));
		# Project name and directory
		group_name_full = group_name+'.'+feat_type+'.'+str(n_clusters)+'cl.'+str(batch_size)+'batch.'+str(n_iters)+'iter.'+str(n_feats)+'feat.'+str(n_timeframes)+'frames';


		##############################
		# CLASSIFY CODE
		##############################
		if CLASSIFY:
			ctlset_path = '../ctl_files/'

			##############################
			# CREATE TRAINING FEATURES
			##############################

			train_feat_read0 = np.zeros((0,n_feats*n_clusters));
			train_feat_read1 = np.zeros((0,n_feats*n_clusters));
			# loop through control files
			ctlset_name = ctl_train_boat;
			with open(ctlset_path+'/'+ctlset_name) as f: 

				# extract file name, and start and end vectors
				for line in f:
					ctl_name = line.rstrip('\n').replace('feat_out','supervect_out');
					ctl_name = ctl_name.replace('spec',feat_type);
					#print('Reading: ' + ctl_name);

					# read data from control file into 2D matrix form
					feat_read = np.loadtxt(ctl_name).flatten();

					# concatenate
					train_feat_read0 = np.vstack((train_feat_read0,feat_read));
			labels0 = np.ones(np.shape(train_feat_read0)[0],dtype=int)*BOAT_N;
			labels0 = labels0.astype(int);

			ctlset_name = ctl_train_heli;
			with open(ctlset_path+'/'+ctlset_name) as f: 

				# extract file name, and start and end vectors
				for line in f:
					ctl_name = line.rstrip('\n').replace('feat_out','supervect_out');
					ctl_name = ctl_name.replace('spec',feat_type);
					#print('Reading: ' + ctl_name);

					# read data from control file into 2D matrix form
					feat_read = np.loadtxt(ctl_name).flatten();

					# concatenate
					train_feat_read1 = np.vstack((train_feat_read1,feat_read));

			labels1 = np.ones(np.shape(train_feat_read1)[0],dtype=int)*HELI_N;
			labels1 = labels1.astype(int);

			train_feat_read = np.vstack((train_feat_read0,train_feat_read1));
			labels_train = np.concatenate((labels0,labels1));


			##############################
			# TRAINING SVM CODE
			##############################
			with open('../trained_svm/' + group_name_full +'.csv','w') as f:
				print('\nTraining SVM\n');
			if svm_type == 'full':
				clf = svm.SVC(C=C,kernel=kernel);
				clf.fit(train_feat_read, labels_train);  
			elif svm_type == 'sgd':
				
				clf = SGDClassifier();
				n_samples = np.shape(train_feat_read)[0];
				n_batches = int(np.floor(1.0*n_samples/batch_size));
				print('	n_batches: '+str(n_batches));
				for n_iter in range(int(n_iters)):
					#print('train iter ' + str(n_iter));
					n_shuffle = range(n_samples);
					random.shuffle(n_shuffle);
					
					if n_batches != 0:
						for n_batch in range(n_batches):
							#print('	batch ' + str(n_batch));
							start_ind = n_batch*batch_size;
							end_ind = start_ind+batch_size-1;
							inds = n_shuffle[start_ind:end_ind];
							clf.partial_fit(train_feat_read[inds,:], labels_train[inds],[BOAT_N,HELI_N]);
						inds = n_shuffle[end_ind:-1];
						clf.partial_fit(train_feat_read[inds,:], labels_train[inds],[BOAT_N,HELI_N]);
					else:
						clf.partial_fit(train_feat_read,labels_train,[BOAT_N,HELI_N]);
				with open('../trained_svm/' + group_name_full +'.csv','w') as f:
					np.savetxt(f,np.asarray(clf.coef_));

				

			##############################
			# TEST AND PRINT ACCURACY
			##############################
			# calculate accuracy over time frames
			print('TRAINING ACC')
			test_feat_read = train_feat_read;
			true_labels = labels_train;
			n_labels = np.size(true_labels);
			predicted_labels = clf.predict(test_feat_read);
			predicted_labels = predicted_labels.astype(int);
			acc_timeframe = 1.*np.sum(predicted_labels==true_labels)/n_labels;
			cf_mat = np.array([[0,0],[0,0]]);
			for n in range(n_labels):
				cf_mat[true_labels[n],predicted_labels[n]] = cf_mat[true_labels[n],predicted_labels[n]]+1;


			with open('../results/'+group_name_full + '.train.frame.metrics','w') as f:
				f.write('accuracy (time frame): ' + str(acc_timeframe)+'\n');
				print('accuracy (time frame): ' + str(acc_timeframe));
				f.write('confusion matrix: '+'\n');
				print('confusion matrix: ');
				print(cf_mat);
				f.write(str(cf_mat[0,:])+'\n');
				f.write(str(cf_mat[1,:])+'\n');
				f.write('true/predicted'+'\n');
				for n in range(n_labels):
					f.write(str(true_labels[n])+str(predicted_labels[n]) + ' ' + '\n');














			##############################
			# CREATE TESTING FEATURES (testing acc)
			##############################
			test_feat_read0 = np.zeros((0,n_feats*n_clusters));
			test_feat_read1 = np.zeros((0,n_feats*n_clusters));
			# loop through control files
			ctlset_name = ctl_test_boat;
			with open(ctlset_path+'/'+ctlset_name) as f: 

				# extract file name, and start and end vectors
				for line in f:
					ctl_name = line.rstrip('\n').replace('feat_out','supervect_out');
					ctl_name = ctl_name.replace('spec',feat_type);
					#print('Reading: ' + ctl_name);

					# read data from control file into 2D matrix form
					feat_read = np.loadtxt(ctl_name).flatten();

					# concatenate
					test_feat_read0 = np.vstack((test_feat_read0,feat_read));
			labels0 = np.ones(np.shape(test_feat_read0)[0],dtype=int)*BOAT_N;
			labels0 = labels0.astype(int);
			ctlset_comb = ctlset_name;
			# loop through control files
			ctlset_name = ctl_test_heli;
			ctlset_comb = ctlset_comb+ctlset_name;
			with open(ctlset_path+'/'+ctlset_name) as f: 

				# extract file name, and start and end vectors
				for line in f:
					ctl_name = line.rstrip('\n').replace('feat_out','supervect_out');
					ctl_name = ctl_name.replace('spec',feat_type);
					#print('Reading: ' + ctl_name);

					# read data from control file into 2D matrix form
					feat_read = np.loadtxt(ctl_name).flatten();

					# concatenate
					test_feat_read1 = np.vstack((test_feat_read1,feat_read));
			labels1 = np.ones(np.shape(test_feat_read1)[0],dtype=int)*HELI_N;
			labels1 = labels1.astype(int);

			##############################
			# TEST AND PRINT ACCURACY
			##############################
			# calculate accuracy over time frames
			print('TESTING ACC')
			test_feat_read = np.vstack((test_feat_read0,test_feat_read1));
			true_labels = np.concatenate((labels0,labels1));
			n_labels = np.size(true_labels);
			predicted_labels = clf.predict(test_feat_read);
			predicted_labels = predicted_labels.astype(int);
			acc_timeframe = 1.*np.sum(predicted_labels==true_labels)/n_labels;
			cf_mat = np.array([[0,0],[0,0]]);
			for n in range(n_labels):
				cf_mat[true_labels[n],predicted_labels[n]] = cf_mat[true_labels[n],predicted_labels[n]]+1;

			with open('../results/'+group_name_full + '.test.frame.metrics','w') as f:
				f.write('accuracy (time frame): ' + str(acc_timeframe)+'\n');
				print('accuracy (time frame): ' + str(acc_timeframe));
				f.write('confusion matrix: '+'\n');
				print('confusion matrix: ');
				print(cf_mat);
				f.write(str(cf_mat[0,:])+'\n');
				f.write(str(cf_mat[1,:])+'\n');
				f.write('true/predicted'+'\n');
				for n in range(n_labels):
					f.write(str(true_labels[n])+str(predicted_labels[n]) + ' ' + '\n');

	