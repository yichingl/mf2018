import numpy as np
from sklearn import linear_model # for logistic regression
import features_helpers as fh # for generating features

# run one epoch and visit every training example to update classifier
def run_epoch(clf, train_filenames, train_labels, classes, sr, clip_duration, mtl_eng):
	n_examples = len(train_filenames);

	# for every random example
	n_shuffle = np.arange(n_examples);
	np.random.shuffle(n_shuffle);
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

	# store data
	clf_dict = {};
	clf_dict["clf"] = clf;
	#clf_dict["train_acc"] = train_acc;

	return clf_dict;

# train_final: train final classifier
def train_final(clf, train_filenames, train_labels, classes, sr, clip_duration, mtl_eng, stop_ne):
	filenames = train_filenames;
	labels = np.append(train_labels);

	n_examples = len(labels);

	for ne in range(stop_ne+1):

		# for every random example
		n_shuffle = range(n_examples);
		random.shuffle(n_shuffle);
		for n in n_shuffle:

			# select current example
			cur_filename = filenames[n];
			cur_label = labels[n];

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

	# store data
	clf_dict = {};
	clf_dict["clf"] = clf;
	#clf_dict["train_acc"] = train_acc;

	return clf_dict;

# fit_data: performs sgd partial fit
def fit_data(clf, features, labels, classes):
	clf.partial_fit(features, labels, classes);
	return clf;

def calc_finalacc(clf, filenames, true_labels_in, sr, clip_duration, mtl_eng):
	n_files = len(filenames);
	n_correct = 0;
	confusion_matrix = np.zeros((2,2));

	for n_file in range(n_files):
		filename = filenames[n_file];
		true_label_in = true_labels_in[n_file];
		data_dict = fh.gen_features_single(filename, true_label_in, sr, clip_duration, mtl_eng);
		features = data_dict["features"];
		true_labels = data_dict["labels"];

		# count number of features
		n_features = np.size(true_labels);
		n_features_half = np.round(n_features/2);

		predicted_labels = classify(clf, features);
		ones_labels_count = np.count_nonzero(predicted_labels);

		# get majority vote for predicted label
		if ones_labels_count>=n_features_half:
			predicted_label = 1;
		else:
			predicted_label = 0;

		# get true label
		true_label = true_labels[0]; #all the same

		# fill in confusion matrix
		confusion_matrix[predicted_label][true_label] += 1;

		print("		true label" + str(true_label) + " predicted as " + str(predicted_label));

		if predicted_label == true_label:
			n_correct += 1;

	print(str(n_correct)+'/'+str(n_files)+' correct');
	print('confusion matrix: ');
	print(confusion_matrix);
	return (1.0*n_correct/n_files);


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
def train_classifier(train_filenames, train_labels, classes, n_epochs, min_epochs, sr, clip_duration, mtl_eng):

	# Initialize classifier
	clf = linear_model.SGDClassifier();

	# Initialize accuracy vectors
	train_acc = np.zeros(n_epochs);
	train_acc_ne = 0.0;

	# For every epoch
	for ne in range(n_epochs):

		print('*** Running epoch ' + str(ne));

		# Run sgd and update classifier
		epoch_dict = run_epoch(clf, train_filenames, train_labels, classes, sr, clip_duration, mtl_eng);
		clf = epoch_dict["clf"];
		#train_acc_ne = epoch_dict["train_acc"];

		#train_acc[ne] = train_acc_ne;

		#print(' ** training accuarcy ' + str(train_acc_ne));

	final_dict = train_final(clf, train_filenames, train_labels, classes, sr, clip_duration, mtl_eng, best_stop_ne);

	clf_out = final_dict["clf"];
	clf_dict = {};
	clf_dict["clf"] = clf_out; 
	#clf_dict["train_acc"] = train_acc;
	return clf_dict;
