import numpy as np;
from sklearn import svm;
#from sklearn.svm import LinearSVC;
#from sklearn.linear_model import SGDClassifier;
import random; # for shuffling
import anurag_cnn.feat_extractor as cnn # to import cnn
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
def gen_cnn(feat, use_gpu=False):
    
    features = cnn.main(feat, use_gpu);
    
    return features;

##############################
# MAIN CODE
##############################

# USER INPUT
# Group name
group_name = 'allcat_bal';
# Define train and test names for each class
boat_train_ctl_name = 'allcat_bal.boat';
heli_train_ctl_name = 'allcat_bal.heli';
boat_test_ctl_name = 'allcat_bal.boat';
heli_test_ctl_name = 'allcat_bal.heli';
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
	# For CNN
	gensupervect_train_ctl_names = [boat_train_ctl_name, heli_train_ctl_name];
	gensupervect_test_ctl_names = [boat_test_ctl_name, heli_test_ctl_name];

# Reading spectrogram path
feat_out = '/home/yichingl/mf/feat_out/';

# Writing CNN_feat path
supervect_out = '/home/yichingl/mf/CNN_feat_out/';

feat_types = ['mel']; #always mel for this CNN, which was trained using melspectrograms


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
	  
	elif feat_type == 'mel':
		#n_feats = 128;
		#n_timeframes = 173;

		n_feats = 128;
		n_timeframes = 173;
	elif feat_type == 'mod':
		n_feats = 128;
		n_timeframes = 320;

	####################################
	# LOAD CNN CODE/MAKE CNN_FEAT CODE
	####################################

	if MAKE_CNN_FEATS:

		# generate CNN features from all data
		# for ctl_group_name_supervect in gensupervect_train_ctl_names:
		# 	ctl_path = '../ctl_files/' + ctl_group_name_supervect +'.train.ctl';
		# 	# loop through control file
		# 	with open(ctl_path) as f: 
		# 		# extract file name, and start and end vectors
		# 		for line in f:
		# 			ctl_name = line.rstrip('\n');
		# 			ctl_name = ctl_name.replace('spec',feat_type);
		# 			print('Reading (gen CNN feat-train): ' + ctl_name);

		# 			# read data from control file into 2D matrix form
		# 			feat_read = np.loadtxt(ctl_name);

		# 			# generate CNN feature from specgram feature
		# 			cnn_feat = gen_cnn(feat_read);

		# 			# save means
		# 			supervect_out_path = ctl_name.replace('feat_out','CNN_feat_out');
		# 			with open(supervect_out_path, 'w') as f:
		# 				np.savetxt(f,np.asarray(cnn_feat));

		for ctl_group_name_supervect in gensupervect_test_ctl_names: # boat, heli
			ctl_path = '../ctl_files/' + ctl_group_name_supervect +'.test.ctl';
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

	if LEARN_PARAMS:

		xval_num_splits = 3; #num folds?

		param_grid = [
		  {'C': np.logspace(-8.,12,num=11,base=2.0), 'gamma': ['auto'], 'kernel': ['linear']},
		  {'C': np.logspace(-8.,12,num=11,base=2.0), 'gamma': np.append(np.logspace(-8.,3,num=7,base=2.0),[1./n_feats]), 'kernel': ['rbf']},
		 ]

		if TEST:
			param_grid = [
			  {'C': [1,100], 'gamma': ['auto'], 'kernel': ['linear']},
			  {'C': [1,100], 'gamma': [0.1,'auto'], 'kernel': ['rbf']},
			 ]

		# for x-val
		split_obj = GroupShuffleSplit(test_size=.33, n_splits=xval_num_splits);

		# Project name and directory
		group_name_full = group_name+'.'+feat_type+'.'+'batch.'+str(n_feats)+'feat.'+str(n_timeframes)+'frames';
		best_params_filepath = '../new_gmm_results/' + group_name_full +'.best_params.txt';
		# Pretty print results filepath
		pretty_results_filepath = '../new_gmm_results/' + group_name_full +'.results.txt';

		with open('../new_cnn_results/' + group_name_full +'.results.txt','w') as f:
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
		xval_num_splits = 5.;

		print('now searching for best params');
		best_avg_mcc = 0;
		best_C = 0;
		best_kernel = None;
		best_gamma = 0;
				
		with open('../new_cnn_results/' + group_name_full +'.results.txt','a') as f:

			# try each classifier parameter set
			for param_type in param_grid:
				kernel_type = param_type['kernel'][0];
				print('kernel type: ' + kernel_type +'\n');
				f.write('kernel type: ' + kernel_type +'\n');

				# for each C parameter 
				for C in param_type['C']:
					f.write('	C: '+str(C)+'\n');
					print('	C: '+str(C)+'\n');

					for gamma in param_type['gamma']:

						f.write('		gamma: '+str(gamma)+'\n');
						print('		gamma: '+str(gamma)+'\n');

						# reset accuracy count for each method
						avg_accuracy = 0;
						avg_mcc = 0;

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

							# for each subset train SVM and make predictions
							clf = svm.SVC(C=C,kernel=kernel_type,gamma=gamma);
							clf.fit(train_feat_subset, train_labels_subset);  
							predicted_test_labels = clf.predict(test_feat_subset);

							# compute performance metrics
							precision = svm_h.calc_precision(test_labels_subset,predicted_test_labels);
							recall = svm_h.calc_recall(test_labels_subset,predicted_test_labels);
							f1 = svm_h.calc_f1(test_labels_subset,predicted_test_labels);
							mcc = svm_h.calc_mcc(test_labels_subset,predicted_test_labels);
							accuracy = svm_h.calc_accuracy(test_labels_subset,predicted_test_labels);
							cf_mat = svm_h.calc_cfmat(test_labels_subset,predicted_test_labels);

							f.write('			Trial MCC: '+ str(mcc)+'\n');
							f.write('				accuracy: ' + str(accuracy)+'\n');
							f.write('				precision: ' + str(precision)+'\n');
							f.write('				recall: ' + str(recall)+'\n');
							f.write('				f1: ' + str(f1)+'\n');
							f.write('				confusion matrix: '+'\n');
							f.write('					'+str(cf_mat[0,:])+'\n');
							f.write('					'+str(cf_mat[1,:])+'\n');
							print('			Trial MCC: '+ str(mcc)+'\n');
							print('				accuracy: ' + str(accuracy)+'\n');
							print('				precision: ' + str(precision)+'\n');
							print('				recall: ' + str(recall)+'\n');
							print('				f1: ' + str(f1)+'\n');
							print('				confusion matrix: '+'\n');
							print('					'+str(cf_mat[0,:])+'\n');
							print('					'+str(cf_mat[1,:])+'\n');

							# update accuracy and sum
							avg_accuracy = avg_accuracy + accuracy
							avg_mcc = avg_mcc + mcc;

						# calcuate average accuracy
						avg_accuracy = 1.*avg_accuracy/xval_num_splits;
						f.write('			Average Accuracy: ' + str(avg_accuracy)+'\n');
						print('			Average Accuracy: ' + str(avg_accuracy)+'\n');

						# calcualate average mcc
						avg_mcc = 1.*avg_mcc/xval_num_splits;
						f.write('			Average MCC: ' + str(avg_mcc)+'\n\n');
						print('			Average MCC: ' + str(avg_mcc)+'\n\n');

						# update best accuracy
						if avg_mcc > best_avg_mcc:
							best_avg_mcc = avg_mcc;
							best_avg_accuracy = avg_accuracy;
							best_C = C;
							best_gamma = gamma
							best_kernel = kernel_type;

		print('best kernel: ' + best_kernel+'\n');
		print('best C: ' + str(best_C)+'\n');
		print('best gamma: ' + str(best_gamma)+'\n');
		print('best avg accuracy: ' + str(best_avg_accuracy)+'\n');
		print('best avg mcc: ' + str(best_avg_mcc)+'\n')
		with open(pretty_results_filepath,'a') as f:
			f.write('\n\n');
			f.write('==================================\n')
			f.write('best kernel: ' + best_kernel+'\n');
			f.write('best C: ' + str(best_C)+'\n');
			f.write('best gamma: ' + str(best_gamma)+'\n');
			f.write('best avg accuracy: ' + str(best_avg_accuracy)+'\n');
			f.write('best avg mcc: ' + str(best_avg_mcc)+'\n');
		

		print('Writing best params to: '+best_params_filepath);
		with open(best_params_filepath,'w') as f:
			f.write(best_kernel+'\n');
			f.write(str(best_C)+'\n')
			f.write(str(best_gamma)+'\n');




	if TRAIN_SVM:

		# preprocess data (train and test sets)
		scaler = preprocessing.StandardScaler().fit(train_feat_read);
		train_feat_read = scaler.transform(train_feat_read);

		# general svm
		clf = svm.SVC(C=best_C,kernel=best_kernel, gamma=best_gamma);
		clf.fit(train_feat_read, labels_train); 

		# make predictions
		predicted_train_labels = clf.predict(train_feat_read);

		# calculate performance metrics
		precision = svm_h.calc_precision(labels_train,predicted_train_labels);
		recall = svm_h.calc_recall(labels_train,predicted_train_labels);
		f1 = svm_h.calc_f1(labels_train,predicted_train_labels);
		mcc = svm_h.calc_mcc(labels_train,predicted_train_labels);
		accuracy = svm_h.calc_accuracy(labels_train,predicted_train_labels);
		cf_mat = svm_h.calc_cfmat(labels_train,predicted_train_labels);

		print('final training accuracy: ' + str(accuracy)+'\n');
		print('final training mcc: ' + str(mcc)+'\n');
		with open('../new_gmm_results/' + group_name_full +'.results.txt','a') as f:
			f.write('\n\n');
			f.write('==================================\n');
			f.write('Training Results');
			f.write('	MCC: '+ str(mcc)+'\n');
			f.write('	accuracy: ' + str(accuracy)+'\n');
			f.write('	precision: ' + str(precision)+'\n');
			f.write('	recall: ' + str(recall)+'\n');
			f.write('	f1: ' + str(f1)+'\n');
			f.write('	confusion matrix: '+'\n');
			f.write('		'+str(cf_mat[0,:])+'\n');
			f.write('		'+str(cf_mat[1,:])+'\n');
		np.savetxt('../new_cnn_results/' + group_name_full +'.train_predicted.txt',predicted_train_labels);

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

		# predict labels
		predicted_test_labels = clf.predict(test_feat_read);

		# calculate performance metrics
		precision = svm_h.calc_precision(labels_test,predicted_test_labels);
		recall = svm_h.calc_recall(labels_test,predicted_test_labels);
		f1 = svm_h.calc_f1(labels_test,predicted_test_labels);
		mcc = svm_h.calc_mcc(labels_test,predicted_test_labels);
		accuracy = svm_h.calc_accuracy(labels_test,predicted_test_labels);
		cf_mat = svm_h.calc_cfmat(labels_test,predicted_test_labels);

		print('final testing accuracy: ' + str(accuracy)+'\n');
		print('final testing mcc: ' + str(mcc)+'\n');
		with open('../new_gmm_results/' + full_name +'.results.txt','a') as f:
			f.write('\n\n');
			f.write('==================================\n')
			f.write('Test Results');
			f.write('	MCC: '+ str(mcc)+'\n');
			f.write('	accuracy: ' + str(accuracy)+'\n');
			f.write('	precision: ' + str(precision)+'\n');
			f.write('	recall: ' + str(recall)+'\n');
			f.write('	f1: ' + str(f1)+'\n');
			f.write('	confusion matrix: '+'\n');
			f.write('		'+str(cf_mat[0,:])+'\n');
			f.write('		'+str(cf_mat[1,:])+'\n');
		np.savetxt('../new_cnn_results/' + group_name_full +'.test_predicted.txt',predicted_test_labels);
