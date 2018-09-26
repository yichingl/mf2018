import numpy as np;
#from sklearn import svm;
from sklearn.svm import LinearSVC;
#from sklearn.linear_model import SGDClassifier;
import random; # for shuffling
#import anurag_cnn.feat_extractor as cnn # to import cnn
#from anurag_cnn.feat_extractor import usegpu
from sklearn import preprocessing
import svm_helper as svm_h
from sklearn.model_selection import GroupShuffleSplit

TEST = 0;

BOAT_N = 0;
HELI_N = 1;
##############################
# TRAIN CNN HELPER CODE
##############################
# gen_cnn: use a pre-trained CNN to extract features
def gen_cnn(feat, use_gpu=True):
    
    features = cnn.main(feat, use_gpu);
    
    return features;

##############################
# MAIN CODE
##############################

# USER INPUT
# Group name
group_name = 'sel';
boat_train_ctl_name = 'sel.boatA';
heli_train_ctl_name = 'sel.heliA';
boat_test_ctl_name = 'sel.boatA';
heli_test_ctl_name = 'sel.heliA';
# For generating CNN_feat
gensupervect_train_ctl_names = [boat_train_ctl_name, heli_train_ctl_name];
gensupervect_test_ctl_names = [boat_test_ctl_name, heli_test_ctl_name];

if TEST:
	print('TESTING');
	group_name = 'test';
	# For SVM
	boat_train_ctl_name = 'test.boat';
	heli_train_ctl_name = 'test.heli';
	boat_test_ctl_name = 'test.boat';
	heli_test_ctl_name = 'test.heli';
	# For GMM
	gensupervect_train_ctl_names = [boat_train_ctl_name, heli_train_ctl_name];
	gensupervect_test_ctl_names = [boat_test_ctl_name, heli_test_ctl_name];

# Reading spectrogram path
feat_out = '/home/yichingl/mf/feat_out/';

# Writing CNN_feat path
supervect_out = '/home/yichingl/mf/CNN_feat_out/';

feat_type = 'mel'; #always mel for this CNN, which was trained using melspectrograms

# Define tasks
MAKE_CNN_FEATS = 0;
#############
LEARN_PARAMS = 1;
TRAIN_SVM = 1; # must use learn_params together? for now
#############
CLASSIFY = 1;

# Define audio params
sr = 44100; # sampling rate
clip_duration = 2; # audio clips will split into segments of this duration (in seconds)


# Define CNN params
use_gpu = True;

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

	n_feats = 133;
	n_timeframes = 669;

	#n_feats = 1025;
	#n_timeframes = 173;
  
elif feat_type == 'mel':
	#n_feats = 128;
	#n_timeframes = 173;

	n_feats = 133;
	n_timeframes = 669;
elif feat_type == 'mod':
	n_feats = 128;
	n_timeframes = 320;

####################################
# LOAD CNN CODE/MAKE CNN_FEAT CODE
####################################

if MAKE_CNN_FEATS:

	# generate CNN features from all data
	for ctl_group_name_supervect in gensupervect_train_ctl_names:
		ctl_path = '../ctl_files/' + ctl_group_name_supervect +'.train.ctl';
		# loop through control file
		with open(ctl_path) as f: 
			# extract file name, and start and end vectors
			for line in f:
				ctl_name = line.rstrip('\n');
				ctl_name = ctl_name.replace('spec',feat_type);
				print('Reading (gen CNN feat-train): ' + ctl_name);

				# read data from control file into 2D matrix form
				feat_read = np.loadtxt(ctl_name);

				# generate CNN feature from specgram feature
				cnn_feat = gen_cnn(feat_read);

				# save means
				supervect_out_path = ctl_name.replace('feat_out','CNN_feat_out');
				with open(supervect_out_path, 'w') as f:
					np.savetxt(f,np.asarray(cnn_feat));

	ctl_path = '../ctl_files/' + group_name +'.test.ctl';
	# loop through control file
	with open(ctl_path) as f: 
		# extract file name, and start and end vectors
		for line in f:
			ctl_name = line.rstrip('\n');
			ctl_name = ctl_name.replace('spec',feat_type);
			print('Reading (gen CNN feat-test): ' + ctl_name);

			# read data from control file into 2D matrix form
			feat_read = np.loadtxt(ctl_name);

			# generate CNN feature from specgram feature
			cnn_feat = gen_cnn(feat_read);

			# save means
			supervect_out_path = ctl_name.replace('feat_out','CNN_feat_out');
			with open(supervect_out_path, 'w') as f:
				np.savetxt(f,np.asarray(cnn_feat));

##############################
# LEARN PARAMETERS CODE
##############################

# define grid to explore using x-validation
if LEARN_PARAMS:

	# param_grid = [
	#   {'C': [1, 10, 100, 1000], 'kernel': ['linear']},
	#   {'C': [1, 10, 100, 1000], 'gamma': [0.001, 0.0001], 'kernel': ['rbf']},
	#  ]

	# param_grid = [
	#   {'C': [50, 75, 100, 125, 150, 250, 500, 750], 'kernel': ['linear']},
	#  ]
	param_grid = [
	  {'C': np.logspace(-3, 2, 6), 'kernel': ['linear']},
	 ]

	# for x-val
	split_obj = GroupShuffleSplit(test_size=.33, n_splits=3);

	for batch_size in batch_size_opts:
		for n_iters in n_iters_opts:
			print('\n\n\nbatch_size: '+str(batch_size)+', n_iters: '+str(n_iters));
			# Project name and directory
			group_name_full = group_name+'.'+feat_type+'.'+'batch.'+str(n_iters)+'iter.'+str(n_feats)+'feat.'+str(n_timeframes)+'frames';
		
			
			with open('../trained_svm_cnn/' + group_name_full +'.xval_results.txt','w') as f:
				ctlset_path = '../ctl_files/';

			##############################
			# CREATE TRAINING FEATURES
			##############################

			# read control file for boat training data, create label vector, read group data
			filename = ctlset_path+boat_train_ctl_name;
			train_feat_read_boat = svm_h.read_cnn_features(filename+'.train.ctl',feat_type);
			labels_boat = svm_h.create_labels(train_feat_read_boat, BOAT_N);
			subgroups_boat = svm_h.read_subgroups(filename+'.train_subgroups.ctl');

			# read control file for heli training data, create label vector, read group data
			filename = ctlset_path+heli_train_ctl_name;
			train_feat_read_heli = svm_h.read_cnn_features(filename+'.train.ctl',feat_type);
			labels_heli = svm_h.create_labels(train_feat_read_heli, HELI_N);
			subgroups_heli = svm_h.read_subgroups(filename+'.train_subgroups.ctl');

			# concatenate boat and heli data
			train_feat_read = np.vstack((train_feat_read_boat,train_feat_read_heli));
			labels_train = np.concatenate((labels_boat,labels_heli));
			subgroups = subgroups_boat+subgroups_heli;

			print('completed reading data');


			##############################
			# TRAINING SVM CODE
			##############################
			print('now searching for best params');
			best_avg_accuracy = 0;
			best_C = 0;
			best_kernel = None;
			best_cfmat = np.array([[0,0],[0,0]]);

			# try each classifier parameter set
			for param_type in param_grid:
				print('kernel type: ' + param_type['kernel'][0]+'\n');

				# for each C parameter 
				for C in param_type['C']:
					with open('../trained_svm_cnn/' + group_name_full +'.xval_results.txt','a') as f:
						f.write('	C: '+str(C)+'\n');

						# reset accuracy count for each method
						avg_accuracy = 0;
						avg_cf_mat = np.array([[0,0],[0,0]]);

						# generate a new set of data - split into data groups
						for train_ind,test_ind in split_obj.split(X=train_feat_read,y=labels_train,groups=subgroups):
							
							# extract train set
							train_idx = np.asarray(train_ind,dtype=int)
							train_feat_subset = train_feat_read[train_idx,:];
							train_labels_subset = labels_train[train_idx];
							# extract test set
							test_idx = np.asarray(test_ind,dtype=int)
							test_feat_subset = train_feat_read[test_idx,:];
							test_labels_subset = labels_train[test_idx];

							# preprocess data (train and test sets)
							scaler = preprocessing.StandardScaler().fit(train_feat_subset);
							train_feat_subset = scaler.transform(train_feat_subset);
							test_feat_subset = scaler.transform(test_feat_subset);

							# for each subset train SVM and calculate accuracy
							# general svm
							#clf = svm.SVC(C=C,kernel=param_type['kernel'][0]);
							#clf.fit(train_feat_subset, train_labels_subset);  
							# linear svm
							clf = LinearSVC(C=C);
							clf.fit(train_feat_subset, train_labels_subset);  

							predicted_test_labels = clf.predict(test_feat_subset);
							accuracy = svm_h.calc_accuracy(test_labels_subset,predicted_test_labels);
							f.write('		accuracy: ' + str(accuracy)+'\n');

							# calculate and update confusion matrix
							cf_mat = svm_h.calc_cfmat(test_labels_subset,predicted_test_labels);

							# update accuracy sum
							avg_accuracy = avg_accuracy + accuracy

							# update confusion matrix sum
							avg_cf_mat = avg_cf_mat + cf_mat

						# calcuate average accuracy
						avg_accuracy = 1.*avg_accuracy/3.;
						f.write('		average accuracy: ' + str(avg_accuracy)+'\n');

						# calcuate average confusion matrix
						avg_cf_mat = np.multiply(1./3, avg_cf_mat);
						f.write('	avg confusion matrix: '+'\n');
						f.write('	'+str(cf_mat[0,:])+'\n');
						f.write('	'+str(cf_mat[1,:])+'\n\n');


						# update best accuracy
						if avg_accuracy > best_avg_accuracy:
							best_avg_accuracy = avg_accuracy;
							best_C = C
							best_kernel = param_type['kernel'][0];
							best_cfmat = avg_cf_mat;
			with open('../trained_svm_cnn/' + group_name_full +'.xval_results.txt','a') as f:
				f.write('\n\n');
				f.write('==================================\n')
				f.write('best kernel: ' + best_kernel+'\n');
				f.write('best C: ' + str(best_C)+'\n');
				f.write('best avg accuracy: ' + str(best_avg_accuracy)+'\n');




if TRAIN_SVM:

	# preprocess data (train and test sets)
	scaler = preprocessing.StandardScaler().fit(train_feat_read);
	train_feat_read = scaler.transform(train_feat_read);

	# for each subset train SVM and calculate accuracy
	

	# general svm
	#clf = svm.SVC(C=best_C,kernel=best_kernel);
	# clf.fit(train_feat_read, labels_train);  
	# linear svm
	clf = LinearSVC(C=C);
	clf.fit(train_feat_read, labels_train); 
	predicted_train_labels = clf.predict(train_feat_read);
	accuracy = svm_h.calc_accuracy(labels_train,predicted_train_labels);
	print('final training accuracy: ' + str(accuracy)+'\n');
	with open('../trained_svm_cnn/' + group_name_full +'.xval_results.txt','a') as f:
		f.write('\n\n');
		f.write('==================================\n')
		f.write('final training accuracy: ' + str(accuracy));

if CLASSIFY:
	# read control file for boat testing data, create label vector, read group data
	filename = ctlset_path+boat_test_ctl_name;
	test_feat_read_boat = svm_h.read_cnn_features(filename+'.test.ctl',feat_type);
	labels_boat = svm_h.create_labels(test_feat_read_boat, BOAT_N);

	# read control file for heli testing data, create label vector, read group data
	filename = ctlset_path+heli_test_ctl_name;
	test_feat_read_heli = svm_h.read_cnn_features(filename+'.test.ctl',feat_type);
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
	with open('../trained_svm_cnn/' + group_name_full +'.xval_results.txt','a') as f:
		f.write('\n\n');
		f.write('==================================\n')
		f.write('final testing accuracy: ' + str(accuracy));
		f.write('\nconfusion matrix: '+'\n');
		f.write(str(cf_mat[0,:])+'\n');
		f.write(str(cf_mat[1,:])+'\n');



			####
			# with open('../trained_svm_cnn/' + group_name_full +'.csv','w') as f:
			# 	print('\nTraining SVM\n');
			# if svm_type == 'full':
			# 	clf = svm.SVC(C=C,kernel=kernel);
			# 	clf.fit(train_feat_read, labels_train);  
			# elif svm_type == 'sgd':
				
			# 	clf = SGDClassifier();
			# 	n_samples = np.shape(train_feat_read)[0];
			# 	n_batches = int(np.floor(1.0*n_samples/batch_size));
			# 	print('	n_batches: '+str(n_batches));
			# 	for n_iter in range(int(n_iters)):
			# 		#print('train iter ' + str(n_iter));
			# 		n_shuffle = range(n_samples);
			# 		random.shuffle(n_shuffle);
					
			# 		if n_batches != 0:
			# 			for n_batch in range(n_batches):
			# 				#print('	batch ' + str(n_batch));
			# 				start_ind = n_batch*batch_size;
			# 				end_ind = start_ind+batch_size-1;
			# 				inds = n_shuffle[start_ind:end_ind];
			# 				clf.partial_fit(train_feat_read[inds,:], labels_train[inds],[BOAT_N,HELI_N]);
			# 			inds = n_shuffle[end_ind:-1];
			# 			clf.partial_fit(train_feat_read[inds,:], labels_train[inds],[BOAT_N,HELI_N]);
			# 		else:
			# 			clf.partial_fit(train_feat_read,labels_train,[BOAT_N,HELI_N]);
			# with open('../trained_svm_cnn/' + group_name_full +'.csv','w') as f:
			# 	np.savetxt(f,np.asarray(clf.coef_));

				

			# ##############################
			# # TEST AND PRINT ACCURACY
			# ##############################
			# # calculate accuracy over time frames
			# print('TRAINING ACC')
			# test_feat_read = train_feat_read;
			# true_labels = labels_train;
			# n_labels = np.size(true_labels);
			# predicted_labels = clf.predict(test_feat_read);
			# predicted_labels = predicted_labels.astype(int);
			# acc_timeframe = 1.*np.sum(predicted_labels==true_labels)/n_labels;
			# cf_mat = np.array([[0,0],[0,0]]);
			# for n in range(n_labels):
			# 	cf_mat[true_labels[n],predicted_labels[n]] = cf_mat[true_labels[n],predicted_labels[n]]+1;


			# with open('../cnn_results/'+group_name_full + '.train.frame.metrics','w') as f:
			# 	f.write('accuracy (time frame): ' + str(acc_timeframe)+'\n');
			# 	print('accuracy (time frame): ' + str(acc_timeframe));
			# 	f.write('confusion matrix: '+'\n');
			# 	print('confusion matrix: ');
			# 	print(cf_mat);
			# 	f.write(str(cf_mat[0,:])+'\n');
			# 	f.write(str(cf_mat[1,:])+'\n');
			# 	f.write('true/predicted'+'\n');
			# 	for n in range(n_labels):
			# 		f.write(str(true_labels[n])+str(predicted_labels[n]) + ' ' + '\n');














			# ##############################
			# # CREATE TESTING FEATURES (testing acc)
			# ##############################
			# test_feat_read0 = np.zeros((0,1024));
			# test_feat_read1 = np.zeros((0,1024));
			# # loop through control files
			# ctlset_name = boat_test_ctl_name+'.test.ctl';
			# with open(ctlset_path+'/'+ctlset_name) as f: 

			# 	# extract file name, and start and end vectors
			# 	for line in f:
			# 		ctl_name = line.rstrip('\n').replace('feat_out','CNN_feat_out');
			# 		ctl_name = ctl_name.replace('spec',feat_type);
			# 		#print('Reading: ' + ctl_name);

			# 		# read data from control file into 2D matrix form
			# 		feat_read = np.loadtxt(ctl_name).flatten();

			# 		# concatenate
			# 		test_feat_read0 = np.vstack((test_feat_read0,feat_read));
			# labels0 = np.ones(np.shape(test_feat_read0)[0],dtype=int)*BOAT_N;
			# labels0 = labels0.astype(int);
			# # loop through control files
			# ctlset_name = heli_test_ctl_name+'.test.ctl';
			# with open(ctlset_path+'/'+ctlset_name) as f: 

			# 	# extract file name, and start and end vectors
			# 	for line in f:
			# 		ctl_name = line.rstrip('\n').replace('feat_out','CNN_feat_out');
			# 		ctl_name = ctl_name.replace('spec',feat_type);
			# 		#print('Reading: ' + ctl_name);

			# 		# read data from control file into 2D matrix form
			# 		feat_read = np.loadtxt(ctl_name).flatten();

			# 		# concatenate
			# 		test_feat_read1 = np.vstack((test_feat_read1,feat_read));
			# labels1 = np.ones(np.shape(test_feat_read1)[0],dtype=int)*HELI_N;
			# labels1 = labels1.astype(int);

			# ##############################
			# # TEST AND PRINT ACCURACY
			# ##############################
			# # calculate accuracy over time frames
			# print('TESTING ACC')
			# test_feat_read = np.vstack((test_feat_read0,test_feat_read1));
			# true_labels = np.concatenate((labels0,labels1));
			# n_labels = np.size(true_labels);
			# predicted_labels = clf.predict(test_feat_read);
			# predicted_labels = predicted_labels.astype(int);
			# acc_timeframe = 1.*np.sum(predicted_labels==true_labels)/n_labels;
			# cf_mat = np.array([[0,0],[0,0]]);
			# for n in range(n_labels):
			# 	cf_mat[true_labels[n],predicted_labels[n]] = cf_mat[true_labels[n],predicted_labels[n]]+1;

			# with open('../cnn_results/'+group_name_full + '.test.frame.metrics','w') as f:
			# 	f.write('accuracy (time frame): ' + str(acc_timeframe)+'\n');
			# 	print('accuracy (time frame): ' + str(acc_timeframe));
			# 	f.write('confusion matrix: '+'\n');
			# 	print('confusion matrix: ');
			# 	print(cf_mat);
			# 	f.write(str(cf_mat[0,:])+'\n');
			# 	f.write(str(cf_mat[1,:])+'\n');
			# 	f.write('true/predicted'+'\n');
			# 	for n in range(n_labels):
			# 		f.write(str(true_labels[n])+str(predicted_labels[n]) + ' ' + '\n');

	