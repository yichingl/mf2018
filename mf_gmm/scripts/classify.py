import numpy as np;
from sklearn import svm;

ctlset_path = '../etc';
ctlset_names_test = ['Xtestclass.boat.ctl','Xtestclass.heli.ctl'];
ctlset_names_train = ['Xtrainclass.boat.ctl','Xtrainclass.heli.ctl'];

n_feats = 104;

train_feat_read0 = np.zeros((0,n_feats));
train_feat_read1 = np.zeros((0,n_feats));
# loop through control files
ctlset_name = 'trainclass.boat.ctl';
with open(ctlset_path+'/'+ctlset_name) as f: 

	# extract file name, and start and end vectors
	for line in f:
		line = line.rstrip('\n').replace('feat_out','supervect_out');
		txt = line.split(" ");
		ctl_name = txt[0];
		start_n = int(txt[1]);
		end_n = int(txt[2]);
		print('Reading: ' + ctl_name);

		# read data from control file into 2D matrix form
		feat_read = np.loadtxt(ctl_name);

		# concatenate
		train_feat_read0 = np.vstack((train_feat_read0,feat_read));
labels0 = np.zeros(np.shape(train_feat_read0)[0]);

ctlset_name = 'trainclass.heli.ctl';
with open(ctlset_path+'/'+ctlset_name) as f: 

	# extract file name, and start and end vectors
	for line in f:
		line = line.rstrip('\n').replace('feat_out','supervect_out');
		txt = line.split(" ");
		ctl_name = txt[0];
		start_n = int(txt[1]);
		end_n = int(txt[2]);
		print('Reading: ' + ctl_name);

		# read data from control file into 2D matrix form
		feat_read = np.loadtxt(ctl_name);

		# concatenate
		train_feat_read1 = np.vstack((train_feat_read1,feat_read));

labels1 = np.ones(np.shape(train_feat_read1)[0]);

train_feat_read = np.vstack((train_feat_read0,train_feat_read1));
labels_train = np.concatenate((labels0,labels1));


clf = svm.SVC();
clf.fit(train_feat_read, labels_train);  




test_feat_read0 = np.zeros((0,n_feats));
test_feat_read1 = np.zeros((0,n_feats));
# loop through control files
ctlset_name = 'testclass.boat.ctl';
with open(ctlset_path+'/'+ctlset_name) as f: 

	# extract file name, and start and end vectors
	for line in f:
		ctl_name = line.rstrip('\n').replace('feat_out','supervect_out');
		print('Reading: ' + ctl_name);

		# read data from control file into 2D matrix form
		feat_read = np.loadtxt(ctl_name);

		# concatenate
		test_feat_read0 = np.vstack((test_feat_read0,feat_read));
labels0 = np.zeros(np.shape(train_feat_read0)[0]);
# loop through control files
ctlset_name = 'testclass.heli.ctl';
with open(ctlset_path+'/'+ctlset_name) as f: 

	# extract file name, and start and end vectors
	for line in f:
		ctl_name = line.rstrip('\n').replace('feat_out','supervect_out');
		print('Reading: ' + ctl_name);

		# read data from control file into 2D matrix form
		feat_read = np.loadtxt(ctl_name);

		# concatenate
		test_feat_read1 = np.vstack((test_feat_read1,feat_read));
labels1 = np.ones(np.shape(train_feat_read1)[0]);


test_feat_read = np.vstack((train_feat_read0,train_feat_read1));
labels_test = np.concatenate((labels0,labels1));
predicted_labels = clf.predict(test_feat_read);
print(predicted_labels);
print(labels_test);
print(1.*np.sum(predicted_labels==labels_test)/np.size(labels_test));
