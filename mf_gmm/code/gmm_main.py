import numpy as np;
from sklearn import svm;
import bob.learn.em
from sklearn import preprocessing
import svm_helper as svm_h
import gmm_helper as gmm_h
from sklearn.model_selection import GroupShuffleSplit

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
# Group name
group_name = 'sel';
boat_train_ctl_name = 'sel.boatA';
heli_train_ctl_name = 'sel.heliB';
boat_test_ctl_name = 'sel.boatA';
heli_test_ctl_name = 'sel.heliA';
# For GMM
traingmm_ctl_names = ['sel.boatA','sel.heliB'];
gensupervect_train_ctl_names = [boat_train_ctl_name, heli_train_ctl_name];
gensupervect_test_ctl_names = [boat_test_ctl_name, heli_test_ctl_name];

if TEST:
	print('TESTING');
	group_name = 'test';
	# For SVM
	boat_train_ctl_name = 'test.boat';
	heli_train_ctl_name = 'test.heli';
	boat_test_ctl_name = 'minitest.boatA';
	heli_test_ctl_name = 'minitest.heliA';
	# For GMM
	traingmm_ctl_names = ['test.boat','test.heli'];
	gensupervect_train_ctl_names = [boat_train_ctl_name, heli_train_ctl_name];
	gensupervect_test_ctl_names = [boat_test_ctl_name, heli_test_ctl_name];

# Reading spectrogram path
feat_out = '/home/yichingl/mf/feat_out/';

# Writing supervect path
supervect_out = '/home/yichingl/mf/supervect_out/';

feat_types = ['mod','spec','mel'];

# Define tasks
TRAIN_GMM = 1;
MAKE_SUPERVECTS = 1;  # run this every time new timeframe or ncl
#############
LEARN_PARAMS = 0;
TRAIN_SVM = 1; # must use learn_params together? for now
CLASSIFY = 1;

# Define audio params
sr = 44100; # sampling rate
clip_duration = 2; # audio clips will split into segments of this duration (in seconds)

# Define gmm params
n_clusters = 64;
max_iter = 5;
covariance_type = 'diag';

# Define normalizing method
norm_method = 'minmax'; # 'standard'
if norm_method == 'minmax':
	scaler = preprocessing.MinMaxScaler();
elif norm_method == 'standard':
	scaler = preprocessing.StandardScaler();

for feat_type in feat_types:
	# hardcode (specgram)
	if feat_type == 'spec':
		#n_feats = 442;
		#n_timeframes = 200;

		#n_feats = 104;
		#n_timeframes = 849;

		n_feats = 133;
		n_timeframes = 669;

		#n_feats = 1025;
		#n_timeframes = 173;
		if LEARN_PARAMS == 0:
			best_kernel = 'linear';
			best_C = 4096.0;
			best_gamma: 'auto';

	elif feat_type == 'mel':
		n_feats = 128;
		n_timeframes = 173;

		if LEARN_PARAMS == 0:
			best_kernel = 'rbf';
			best_C = 4096.0;
			best_gamma = 0.001;

		# n_feats = 133;
		# n_timeframes = 669;
	elif feat_type == 'mod':
		n_feats = 128;
		n_timeframes = 320;
		if LEARN_PARAMS == 0:
			best_kernel = 'rbf';
			best_C = 1024.0;
			best_gamma = 'auto';

	##############################
	# TRAIN GMM USING MLE
	##############################

	# Read train ctl file
	if TRAIN_GMM:

		# Initialize matrix for features
		feat_all = np.zeros((0,n_feats));

		# Generate output save paths
		gmm_means_path = '../trained_gmms/'+group_name+'.'+feat_type+'.cl'+str(n_clusters)+'.gmm.means.csv';
		
		# Read all ctl files
		for ctl_group_name in traingmm_ctl_names:

			# Generate input ctl read path
			ctl_path = '../ctl_files/' + ctl_group_name +'.train.ctl';

			# Read ctl file
			read_feat = gmm_h.read_specgram_features(ctl_path, n_feats, feat_type);

			# Add to full matrix
			feat_all = np.vstack((feat_all, read_feat));

		# Normalize data along feature axis
		scaler.fit(feat_all);
		feat_all = scaler.transform(feat_all);

		## Run K-Means clustering to gather means for training GMM
		kmeans_means = gmm_h.get_kmeans_means(feat_all, n_clusters, n_feats);

		## Initialize and train GMM
		gmm_ubm = gmm_h.train_gmm(feat_all, n_clusters, n_feats); #max iter, covariance type?
		gmm_means = gmm_ubm.means;
		# Save GMM parameters
		with open(gmm_means_path, 'w') as f:
			np.savetxt(f,np.asarray(gmm_means));


		# # User input
		# for ctl_group_name in traingmm_ctl_names:
		# 	# don't use start and end for now
		# 	gmm_weights_path = '../trained_gmms/'+ctl_group_name+'.'+feat_type+'.cl'+str(n_clusters)+'.gmm.weights.csv';
		# 	gmm_means_path = '../trained_gmms/'+ctl_group_name+'.'+feat_type+'.cl'+str(n_clusters)+'.gmm.means.csv';
		# 	ctl_path = '../ctl_files/' + ctl_group_name +'.train.ctl';

		# 	# initialize gmm
		# 	gmm = mixture.GaussianMixture(n_components=n_clusters, max_iter=max_iter, covariance_type=covariance_type);

		# 	# hardcode
		# 	feat_all = np.zeros((0,n_feats));

		# 	# loop through control file
		# 	with open(ctl_path) as f: 

		# 		# extract file name, and start and end vectors
		# 		for ctl_name in f:
		# 			ctl_name = ctl_name.rstrip('\n');
		# 			ctl_name = ctl_name.replace('spec',feat_type);
		# 			print('Reading (train gmm): ' + ctl_name);

		# 			# read data
		# 			x = np.transpose(np.loadtxt(ctl_name));

		# 			# normalize data
		# 			x = x/(np.linalg.norm(x, ord=1, axis=1, keepdims=True)+np.finfo(float).eps);
		# 			print(np.shape(x));

		# 			# read data from file and stack
		# 			feat_all = np.vstack((feat_all,x));

		# 	gmm.fit(feat_all);
		# 	# save gmm model
		# 	with open(gmm_weights_path, 'w') as f:
		# 		np.savetxt(f,np.asarray(gmm.weights_));
		# 	with open(gmm_means_path, 'w') as f:
		# 		np.savetxt(f,np.asarray(gmm.means_));

	##############################
	# MAKE SUPERVECTOR CODE
	##############################
	if MAKE_SUPERVECTS:

		# Load in init params (weights, means) for GMM
		gmm_means_path = '../trained_gmms/'+group_name+'.'+feat_type+'.cl'+str(n_clusters)+'.gmm.means.csv';
		init_means = np.loadtxt(gmm_means_path);

		# Process training set (boat, heli)
		for ctl_group_name_supervect in gensupervect_train_ctl_names: # boat, heli
			ctl_path = '../ctl_files/' + ctl_group_name_supervect +'.train.ctl';
			with open(ctl_path) as f: # open boat,heli ctl file

				for line in f: # for each filepath, generate a supervect and save
					ctl_name = gmm_h.read_ctlname(line,feat_type);
					feat_read = gmm_h.read_specgram_feature(ctl_name);
					supervect = gmm_h.create_supervect(feat_read, n_clusters, n_feats, gmm_means);

					# save means
					supervect_out_path = ctl_name.replace('feat_out','supervect_out');
					with open(supervect_out_path, 'w') as f:
						np.savetxt(f,np.asarray(supervect));

		# Process testing set (boat,heli)
		for ctl_group_name_supervect in gensupervect_test_ctl_names: # boat, heli
			ctl_path = '../ctl_files/' + ctl_group_name_supervect +'.test.ctl';
			with open(ctl_path) as f: 
				for line in f: # for each filepath, generate a supervect and save
					ctl_name = gmm_h.read_ctlname(line,feat_type);
					feat_read = gmm_h.read_specgram_feature(ctl_name);
					supervect = gmm_h.create_supervect(feat_read, n_clusters, n_feats, gmm_means);

					# save means
					supervect_out_path = ctl_name.replace('feat_out','supervect_out');
					with open(supervect_out_path, 'w') as f:
						np.savetxt(f,np.asarray(supervect));



		# # Read in GMMs for each category and agglomerate
		# for ctl_group_name_supervect in gensupervect_train_ctl_names:
		# 	# don't use start and end for now
		# 	gmm_weights_path = '../trained_gmms/'+ctl_group_name_supervect+'.'+feat_type+'.cl'+str(n_clusters)+'.gmm.weights.csv';
		# 	gmm_means_path = '../trained_gmms/'+ctl_group_name_supervect+'.'+feat_type+'.cl'+str(n_clusters)+'.gmm.means.csv';

		# 	try:
		# 		init_weights = np.stack((init_weights,np.loadtxt(gmm_weights_path)),axis=0);
		# 		init_means = np.stack((init_means,np.loadtxt(gmm_means_path)),axis=0);
		# 	except NameError:
		# 		# load gmm
		# 		init_weights = np.loadtxt(gmm_weights_path);
		# 		init_means = np.loadtxt(gmm_means_path);

		# # agglomarate and normalize
		# init_weights = np.mean(init_weights,axis=0);
		# init_weights = init_weights/(np.linalg.norm(init_weights, ord=1)+np.finfo(float).eps);

		# init_means = np.mean(init_means,axis=0);
		# init_means = init_means/(np.linalg.norm(init_means, ord=1, axis=1, keepdims=True)+np.finfo(float).eps);


		# # Convert train data to supervector
		# for ctl_group_name_supervect in gensupervect_train_ctl_names:
		# 	ctl_path = '../ctl_files/' + ctl_group_name_supervect +'.train.ctl';
		# 	# loop through control file
		# 	with open(ctl_path) as f: 
		# 		# extract file name, and start and end vectors
		# 		for line in f:
		# 			ctl_name = line.rstrip('\n');
		# 			ctl_name = ctl_name.replace('spec',feat_type);
		# 			print('Reading (supervector-train): ' + ctl_name);

		# 			# read data from control file into 2D matrix form
		# 			feat_read = np.transpose(np.loadtxt(ctl_name));

		# 			# update supervector gmm once
		# 			gmm = mixture.GaussianMixture(n_components=n_clusters, max_iter=1, covariance_type=covariance_type, weights_init=init_weights, means_init=init_means);
		# 			gmm.fit(feat_read);

		# 			# save means
		# 			supervect_out_path = ctl_name.replace('feat_out','supervect_out');
		# 			with open(supervect_out_path, 'w') as f:
		# 				np.savetxt(f,np.asarray(gmm.means_));

		# ctl_path = '../ctl_files/' + group_name +'.test.ctl';
		# # loop through control file
		# with open(ctl_path) as f: 
		# 	# extract file name, and start and end vectors
		# 	for line in f:
		# 		ctl_name = line.rstrip('\n');
		# 		ctl_name = ctl_name.replace('spec',feat_type);
		# 		print('Reading (supervector-test): ' + ctl_name);

		# 		# read data from control file into 2D matrix form
		# 		feat_read = np.transpose(np.loadtxt(ctl_name));
		# 		print(np.shape(feat_read));

		# 		# update supervector gmm once
		# 		gmm = mixture.GaussianMixture(n_components=n_clusters, max_iter=1, covariance_type=covariance_type, weights_init=init_weights, means_init=init_means);
		# 		gmm.fit(feat_read);

		# 		# save means
		# 		supervect_out_path = ctl_name.replace('feat_out','supervect_out');
		# 		with open(supervect_out_path, 'w') as f:
		# 			np.savetxt(f,np.asarray(gmm.means_));

	##############################
	# LEARN PARAMETERS CODE
	##############################
	# Project name and directory
	group_name_full = group_name+'.'+feat_type+'.'+str(n_clusters)+'cl.'+'batch.'+str(n_feats)+'feat.'+str(n_timeframes)+'frames';
	ctlset_path = '../ctl_files/';

	if LEARN_PARAMS or TRAIN_SVM:
		##############################
		# CREATE TRAINING FEATURES
		##############################

		# read control file for boat training data, create label vector, read group data
		filename = ctlset_path+boat_train_ctl_name;
		train_feat_read_boat = svm_h.read_gmm_features(filename+'.train.ctl',feat_type,n_feats*n_clusters);
		labels_boat = svm_h.create_labels(train_feat_read_boat, BOAT_N);
		subgroups_boat = svm_h.read_subgroups(filename+'.train_subgroups.ctl');

		# read control file for heli training data, create label vector, read group data
		filename = ctlset_path+heli_train_ctl_name;
		train_feat_read_heli = svm_h.read_gmm_features(filename+'.train.ctl',feat_type,n_feats*n_clusters);
		labels_heli = svm_h.create_labels(train_feat_read_heli, HELI_N);
		subgroups_heli = svm_h.read_subgroups(filename+'.train_subgroups.ctl');

		# concatenate boat and heli data
		train_feat_read = np.vstack((train_feat_read_boat,train_feat_read_heli));
		labels_train = np.concatenate((labels_boat,labels_heli));
		subgroups = subgroups_boat+subgroups_heli;

		print('completed reading data');

	# define grid to explore using x-validation
	if LEARN_PARAMS:
		num_splits = 5; #num folds?

		param_grid = [
		  {'C': np.logspace(-8.,12,num=11,base=2.0), 'gamma': ['auto'], 'kernel': ['linear']},
		  {'C': np.logspace(-8.,12,num=11,base=2.0), 'gamma': [0.001, 0.0001,'auto'], 'kernel': ['rbf']},
		 ]

		if TEST:
			param_grid = [
			  {'C': [1,100], 'gamma': ['auto'], 'kernel': ['linear']},
			  {'C': [1,100], 'gamma': [0.1,'auto'], 'kernel': ['rbf']},
			 ]

		# for x-val
		split_obj = GroupShuffleSplit(test_size=.33, n_splits=num_splits);

		
		# Clear file
		with open('../new_gmm_results/' + group_name_full +'.results.txt','w') as f:
			print('');

		


		##############################
		# TRAINING SVM CODE
		##############################
		print('now searching for best params');
		best_avg_accuracy = 0;
		best_C = 0;
		best_gamma = 0;
		best_kernel = None;

		with open('../new_gmm_results/' + group_name_full +'.results.txt','a') as f:

			# try each classifier parameter set
			for param_type in param_grid:
				kernel_type = param_type['kernel'][0];
				print('kernel type: ' + kernel_type +'\n');
				f.write('Kernel Type: '+str(kernel_type)+'\n');

				# for each C parameter 
				for C in param_type['C']:
					print('	C: '+str(C)+'\n');
					f.write('	C: '+str(C)+'\n');

					# for each gamma parameter
					for gamma in param_type['gamma']:
						print('		gamma: '+str(gamma)+'\n');
						f.write('		gamma: '+str(gamma)+'\n');

						# reset accuracy count for each method
						avg_accuracy = 0;

						# generate a new set of data - split into data groups
						for train_ind,test_ind in split_obj.split(X=train_feat_read,y=labels_train,groups=subgroups):
							
							# move this?
							# extract train set
							train_idx = np.asarray(train_ind,dtype=int)
							train_feat_subset = train_feat_read[train_idx,:];
							train_labels_subset = labels_train[train_idx];
							# extract test set
							test_idx = np.asarray(test_ind,dtype=int)
							test_feat_subset = train_feat_read[test_idx,:];
							test_labels_subset = labels_train[test_idx];

							# preprocess data (train and test sets)
							scaler.fit(train_feat_subset);
							train_feat_subset = scaler.transform(train_feat_subset);
							test_feat_subset = scaler.transform(test_feat_subset);

							# for each subset train SVM and calculate accuracy
							clf = svm.SVC(C=C,kernel=kernel_type,gamma=gamma);
							clf.fit(train_feat_subset, train_labels_subset);  
							predicted_test_labels = clf.predict(test_feat_subset);
							accuracy = svm_h.calc_accuracy(test_labels_subset,predicted_test_labels);
							f.write('			accuracy: ' + str(accuracy)+'\n');

							# update accuracy sum
							avg_accuracy = avg_accuracy + accuracy
						# calcuate average accuracy
						avg_accuracy = 1.*avg_accuracy/num_splits;
						f.write('			average accuracy: ' + str(avg_accuracy)+'\n\n');
						print('			average accuracy: ' + str(avg_accuracy)+'\n\n');

						# update best accuracy
						if avg_accuracy > best_avg_accuracy:
							best_avg_accuracy = avg_accuracy;
							best_C = C;
							best_gamma = gamma
							best_kernel = kernel_type;

		with open('../new_gmm_results/' + group_name_full +'.results.txt','a') as f:
			f.write('\n\n');
			f.write('==================================\n')
			f.write('best kernel: ' + best_kernel+'\n');
			f.write('best C: ' + str(best_C)+'\n');
			f.write('best gamma: ' + str(best_gamma)+'\n');
			f.write('best avg accuracy: ' + str(best_avg_accuracy)+'\n');
		print('best kernel: ' + best_kernel+'\n');
		print('best C: ' + str(best_C)+'\n');
		print('best gamma: ' + str(best_gamma)+'\n');
		print('best avg accuracy: ' + str(best_avg_accuracy)+'\n');




	if TRAIN_SVM:

		# preprocess data (train and test sets)
		scaler.fit(train_feat_read);
		train_feat_read = scaler.transform(train_feat_read);

		# for each subset train SVM and calculate accuracy
		clf = svm.SVC(C=best_C,gamma=best_gamma,kernel=best_kernel);
		clf.fit(train_feat_read, labels_train);  
		predicted_train_labels = clf.predict(train_feat_read);
		accuracy = svm_h.calc_accuracy(labels_train,predicted_train_labels);
		cf_mat = svm_h.calc_cfmat(labels_train,predicted_train_labels);
		print('final training accuracy: ' + str(accuracy)+'\n');
		with open('../new_gmm_results/' + group_name_full +'.results.txt','a') as f:
			f.write('\n\n');
			f.write('==================================\n')
			f.write('final training accuracy: ' + str(accuracy));
			f.write('\nconfusion matrix: '+'\n');
			f.write(str(cf_mat[0,:])+'\n');
			f.write(str(cf_mat[1,:])+'\n');
		np.savetxt('../new_gmm_results/' + group_name_full +'.train_predicted.txt',predicted_train_labels);
			

	if CLASSIFY:
		# read control file for boat testing data, create label vector, read group data
		filename = ctlset_path+boat_test_ctl_name;
		test_feat_read_boat = svm_h.read_gmm_features(filename+'.test.ctl',feat_type,n_feats*n_clusters);
		labels_boat = svm_h.create_labels(test_feat_read_boat, BOAT_N);

		# read control file for heli testing data, create label vector, read group data
		filename = ctlset_path+heli_test_ctl_name;
		test_feat_read_heli = svm_h.read_gmm_features(filename+'.test.ctl',feat_type,n_feats*n_clusters);
		labels_heli = svm_h.create_labels(test_feat_read_heli, HELI_N);

		# concatenate boat and heli data
		test_feat_read = np.vstack((test_feat_read_boat,test_feat_read_heli));
		labels_test = np.concatenate((labels_boat,labels_heli));

		# preprocess test data
		test_feat_read = scaler.transform(test_feat_read);

		# predict labels and calculate accuracy
		predicted_test_labels = clf.predict(test_feat_read);
		accuracy = svm_h.calc_accuracy(labels_test,predicted_test_labels);
		cf_mat = svm_h.calc_cfmat(labels_test,predicted_test_labels);
		print('final testing accuracy: ' + str(accuracy)+'\n');
		with open('../new_gmm_results/' + group_name_full +'.results.txt','a') as f:
			f.write('\n\n');
			f.write('==================================\n')
			f.write('final testing accuracy: ' + str(accuracy));
			f.write('\nconfusion matrix: '+'\n');
			f.write(str(cf_mat[0,:])+'\n');
			f.write(str(cf_mat[1,:])+'\n');
		np.savetxt('../new_gmm_results/' + group_name_full +'.test_predicted.txt',predicted_test_labels);
	# ##############################
	# # CLASSIFY CODE
	# ##############################
	# for batch_size in batch_size_opts:
	# 	for n_iters in n_iters_opts:
	# 		print('\n\n\nn_clusters: '+str(n_clusters)+', batch_size: '+str(batch_size)+', n_iters: '+str(n_iters));
	# 		# Project name and directory
	# 		group_name_full = group_name+'.'+feat_type+'.'+str(n_clusters)+'cl.'+str(batch_size)+'batch.'+str(n_iters)+'iter.'+str(n_feats)+'feat.'+str(n_timeframes)+'frames';



	# 		if CLASSIFY:
	# 			ctlset_path = '../ctl_files/'

	# 			##############################
	# 			# CREATE TRAINING FEATURES
	# 			##############################

	# 			train_feat_read0 = np.zeros((0,n_feats*n_clusters));
	# 			train_feat_read1 = np.zeros((0,n_feats*n_clusters));
	# 			# loop through control files
	# 			ctlset_name = boat_train_ctl_name+'.train.ctl';
	# 			with open(ctlset_path+'/'+ctlset_name) as f: 

	# 				# extract file name, and start and end vectors
	# 				for line in f:
	# 					ctl_name = line.rstrip('\n').replace('feat_out','supervect_out');
	# 					ctl_name = ctl_name.replace('spec',feat_type);
	# 					# this should actually be gmm type???
	# 					#print('Reading: ' + ctl_name);

	# 					# read data from control file into 2D matrix form
	# 					feat_read = np.loadtxt(ctl_name).flatten();

	# 					# concatenate
	# 					#print(np.shape(train_feat_read0));
	# 					#print(np.shape(feat_read));
	# 					train_feat_read0 = np.vstack((train_feat_read0,feat_read));
	# 			labels0 = np.ones(np.shape(train_feat_read0)[0],dtype=int)*BOAT_N;
	# 			labels0 = labels0.astype(int);

	# 			ctlset_name = heli_train_ctl_name+'.train.ctl';
	# 			with open(ctlset_path+'/'+ctlset_name) as f: 

	# 				# extract file name, and start and end vectors
	# 				for line in f:
	# 					ctl_name = line.rstrip('\n').replace('feat_out','supervect_out');
	# 					ctl_name = ctl_name.replace('spec',feat_type);
	# 					#print('Reading: ' + ctl_name);

	# 					# read data from control file into 2D matrix form
	# 					feat_read = np.loadtxt(ctl_name).flatten();

	# 					# concatenate
	# 					train_feat_read1 = np.vstack((train_feat_read1,feat_read));

	# 			labels1 = np.ones(np.shape(train_feat_read1)[0],dtype=int)*HELI_N;
	# 			labels1 = labels1.astype(int);

	# 			train_feat_read = np.vstack((train_feat_read0,train_feat_read1));
	# 			labels_train = np.concatenate((labels0,labels1));


	# 			##############################
	# 			# TRAINING SVM CODE
	# 			##############################
	# 			with open('../trained_svm_gmm/' + group_name_full +'.csv','w') as f:
	# 				print('\nTraining SVM\n');
	# 			if svm_type == 'full':
	# 				clf = svm.SVC(C=C,kernel=kernel);
	# 				clf.fit(train_feat_read, labels_train);  
	# 			elif svm_type == 'sgd':
					
	# 				clf = SGDClassifier();
	# 				n_samples = np.shape(train_feat_read)[0];
	# 				n_batches = int(np.floor(1.0*n_samples/batch_size));
	# 				print('	n_batches: '+str(n_batches));
	# 				for n_iter in range(int(n_iters)):
	# 					#print('train iter ' + str(n_iter));
	# 					n_shuffle = range(n_samples);
	# 					random.shuffle(n_shuffle);
						
	# 					if n_batches != 0:
	# 						for n_batch in range(n_batches):
	# 							#print('	batch ' + str(n_batch));
	# 							start_ind = n_batch*batch_size;
	# 							end_ind = start_ind+batch_size-1;
	# 							inds = n_shuffle[start_ind:end_ind];
	# 							clf.partial_fit(train_feat_read[inds,:], labels_train[inds],[BOAT_N,HELI_N]);
	# 						inds = n_shuffle[end_ind:-1];
	# 						clf.partial_fit(train_feat_read[inds,:], labels_train[inds],[BOAT_N,HELI_N]);
	# 					else:
	# 						clf.partial_fit(train_feat_read,labels_train,[BOAT_N,HELI_N]);
	# 			with open('../trained_svm_gmm/' + group_name_full +'.csv','w') as f:
	# 				np.savetxt(f,np.asarray(clf.coef_));

					

	# 			##############################
	# 			# TEST AND PRINT ACCURACY
	# 			##############################
	# 			# calculate accuracy over time frames
	# 			print('TRAINING ACC')
	# 			test_feat_read = train_feat_read;
	# 			true_labels = labels_train;
	# 			n_labels = np.size(true_labels);
	# 			predicted_labels = clf.predict(test_feat_read);
	# 			predicted_labels = predicted_labels.astype(int);
	# 			acc_timeframe = 1.*np.sum(predicted_labels==true_labels)/n_labels;
	# 			cf_mat = np.array([[0,0],[0,0]]);
	# 			for n in range(n_labels):
	# 				cf_mat[true_labels[n],predicted_labels[n]] = cf_mat[true_labels[n],predicted_labels[n]]+1;


	# 			with open('../gmm_results/'+group_name_full + '.train.frame.metrics','w') as f:
	# 				f.write('accuracy (time frame): ' + str(acc_timeframe)+'\n');
	# 				print('accuracy (time frame): ' + str(acc_timeframe));
	# 				f.write('confusion matrix: '+'\n');
	# 				print('confusion matrix: ');
	# 				print(cf_mat);
	# 				f.write(str(cf_mat[0,:])+'\n');
	# 				f.write(str(cf_mat[1,:])+'\n');
	# 				f.write('true/predicted'+'\n');
	# 				for n in range(n_labels):
	# 					f.write(str(true_labels[n])+str(predicted_labels[n]) + ' ' + '\n');














	# 			##############################
	# 			# CREATE TESTING FEATURES (testing acc)
	# 			##############################
	# 			test_feat_read0 = np.zeros((0,n_feats*n_clusters));
	# 			test_feat_read1 = np.zeros((0,n_feats*n_clusters));
	# 			# loop through control files
	# 			ctlset_name = boat_test_ctl_name+'.test.ctl';
	# 			with open(ctlset_path+'/'+ctlset_name) as f: 

	# 				# extract file name, and start and end vectors
	# 				for line in f:
	# 					ctl_name = line.rstrip('\n').replace('feat_out','supervect_out');
	# 					ctl_name = ctl_name.replace('spec',feat_type);
	# 					#print('Reading: ' + ctl_name);

	# 					# read data from control file into 2D matrix form
	# 					feat_read = np.loadtxt(ctl_name).flatten();

	# 					# concatenate
	# 					test_feat_read0 = np.vstack((test_feat_read0,feat_read));
	# 			labels0 = np.ones(np.shape(test_feat_read0)[0],dtype=int)*BOAT_N;
	# 			labels0 = labels0.astype(int);
	# 			# loop through control files
	# 			ctlset_name = heli_test_ctl_name+'.test.ctl';
	# 			with open(ctlset_path+'/'+ctlset_name) as f: 

	# 				# extract file name, and start and end vectors
	# 				for line in f:
	# 					ctl_name = line.rstrip('\n').replace('feat_out','supervect_out');
	# 					ctl_name = ctl_name.replace('spec',feat_type);
	# 					#print('Reading: ' + ctl_name);

	# 					# read data from control file into 2D matrix form
	# 					feat_read = np.loadtxt(ctl_name).flatten();

	# 					# concatenate
	# 					test_feat_read1 = np.vstack((test_feat_read1,feat_read));
	# 			labels1 = np.ones(np.shape(test_feat_read1)[0],dtype=int)*HELI_N;
	# 			labels1 = labels1.astype(int);

	# 			print(np.shape(labels0));
	# 			print(np.shape(labels1));

	# 			##############################
	# 			# TEST AND PRINT ACCURACY
	# 			##############################
	# 			# calculate accuracy over time frames
	# 			print('TESTING ACC')
	# 			test_feat_read = np.vstack((test_feat_read0,test_feat_read1));
	# 			true_labels = np.concatenate((labels0,labels1));
	# 			n_labels = np.size(true_labels);
	# 			predicted_labels = clf.predict(test_feat_read);
	# 			predicted_labels = predicted_labels.astype(int);
	# 			acc_timeframe = 1.*np.sum(predicted_labels==true_labels)/n_labels;
	# 			cf_mat = np.array([[0,0],[0,0]]);
	# 			for n in range(n_labels):
	# 				cf_mat[true_labels[n],predicted_labels[n]] = cf_mat[true_labels[n],predicted_labels[n]]+1;

	# 			with open('../gmm_results/'+group_name_full + '.test.frame.metrics','w') as f:
	# 				f.write('accuracy (time frame): ' + str(acc_timeframe)+'\n');
	# 				print('accuracy (time frame): ' + str(acc_timeframe));
	# 				f.write('confusion matrix: '+'\n');
	# 				print('confusion matrix: ');
	# 				print(cf_mat);
	# 				f.write(str(cf_mat[0,:])+'\n');
	# 				f.write(str(cf_mat[1,:])+'\n');
	# 				f.write('true/predicted'+'\n');
	# 				for n in range(n_labels):
	# 					f.write(str(true_labels[n])+str(predicted_labels[n]) + ' ' + '\n');

	