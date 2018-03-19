import numpy as np
from sklearn import linear_model # for logistic regression
import features_helpers as fh # for generating features
import random # for shuffling

# run one epoch and visit every training example to update classifier
def run_epoch(clf, train_filenames, train_labels, valid_filenames, valid_labels, classes, sr, clip_duration, mtl_eng):
	n_examples = len(train_filenames);

	# for every random example
	n_shuffle = range(n_examples);
	random.shuffle(n_shuffle);
	for n in n_shuffle:

		# select current example
		cur_filename = train_filenames[n];
		cur_label = train_labels[n];

		print('   ** running example ' + str(n));
		print('      * filename: ' + cur_filename);
		print('      * label: ' + str(cur_label));

		# generate feature matrix and labels
		cur_data_dict = fh.gen_features_single(cur_filename, cur_label, sr, clip_duration, mtl_eng);
		cur_features = cur_data_dict["features"];
		cur_labels = cur_data_dict["labels"]; 

		# do a partial fit of weights to that feature (SGD)
		clf = fit_data(clf, cur_features, cur_labels, classes);

	# use classifier to get accuracy on training and validation sets
	#print('   ** calculating training accuracy');
	#train_acc = calc_accuracy(clf, train_filenames, train_labels, sr, clip_duration, mtl_eng); # TODO don't continue to use this
	print('   ** calculating valid accuracy');
	valid_acc = calc_accuracy(clf, valid_filenames, valid_labels, sr, clip_duration, mtl_eng);

	# store data
	clf_dict = {};
	clf_dict["clf"] = clf;
	#clf_dict["train_acc"] = train_acc;
	clf_dict["valid_acc"] = valid_acc;

	return clf_dict;

# fit_data: performs sgd partial fit
def fit_data(clf, features, labels, classes):
	clf.partial_fit(features, labels, classes);
	return clf;

# calc_accuracy: calculate accuracy using classifier
def calc_accuracy(clf, filenames, true_labels, sr, clip_duration, mtl_eng):
	# generate feature matrix and labels
	data_dict = fh.gen_features_multiple(filenames, true_labels, sr, clip_duration, mtl_eng);
	features = data_dict["features"];
	true_labels = data_dict["labels"]; 

	mean_accuracy = clf.score(features, true_labels);
	return mean_accuracy;

# classify: classifies data using classifier and returns predicted labels
def classify(clf, features):
	predicted_labels = clf.predict(features);
	return predicted_labels;

# train_classifier: run algorithm for provided number of epochs
def train_classifier(train_filenames, train_labels, valid_filenames, valid_labels, classes, n_epochs, sr, clip_duration, mtl_eng):

	# Initialize classifier
	clf = linear_model.SGDClassifier();

	# Initialize accuracy vectors
	train_acc = np.zeros(n_epochs);
	valid_acc = np.zeros(n_epochs);
	train_acc_ne = 0.0;
	valid_acc_ne = 0.0;

	# For every epoch
	for ne in range(n_epochs):

		print('*** Running epoch ' + str(ne));

		# Stop training if performance is good enough
		if (valid_acc_ne <= 0.63) or (ne <= 5):

			# Run sgd and update classifier
			epoch_dict = run_epoch(clf, train_filenames, train_labels, valid_filenames, valid_labels, classes, sr, clip_duration, mtl_eng);
			clf = epoch_dict["clf"];
			#train_acc_ne = epoch_dict["train_acc"];
			valid_acc_ne = epoch_dict["valid_acc"];

			#train_acc[ne] = train_acc_ne;
			valid_acc[ne] = valid_acc_ne;

			#print(' ** training accuarcy ' + str(train_acc_ne));
			print(' ** validation accuarcy ' + str(valid_acc_ne));

		else:
			#print('*** Stopping training at early epoch ' + str(ne) + ' (train accuracy = ' + str(train_acc_ne) + ', valid accuracy = ' + str(valid_acc_ne) + ')');
			print('*** Stopping training at early epoch ' + str(ne) + ' (valid accuracy = ' + str(valid_acc_ne) + ')');

	clf_dict = {};
	clf_dict["clf"] = clf; 
	#clf_dict["train_acc"] = train_acc;
	clf_dict["valid_acc"] = valid_acc;
	return clf_dict;
