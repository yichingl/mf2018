import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

labels_train = np.loadtxt("labels.csv",dtype='float');
print(np.shape(labels_train));
print(labels_train[2400:2405]);

GMM=1;
# 3d
if GMM:
	feat_embedded = np.loadtxt("gmm-3d-feat_embedded.csv",dtype='float',delimiter=",");
	plt.scatter(feat_embedded[:, 0], feat_embedded[:, 1],  feat_embedded[:, 2],c = labels_train);
	
	fig = plt.figure()
	ax = fig.add_subplot(111,projection='3d');

	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')
	plt.show();
	plt.savefig('../vis/tsne-gmm-spec-3d.png');
	plt.close()

	plt.scatter(feat_embedded[:, 2], feat_embedded[:, 0], c = labels_train);
	plt.tight_layout();
	plt.show();
	plt.savefig('../vis/tsne-gmm-spec-3d-zx.png');
	plt.close();

	plt.scatter(feat_embedded[:, 2], feat_embedded[:, 1], c = labels_train);
	plt.tight_layout();
	plt.show();
	plt.savefig('../vis/tsne-gmm-spec-3d-zy.png');
	plt.close();
else:
	feat_embedded = np.loadtxt("cnn-3d-feat_embedded.csv",dtype='float',delimiter=",");
	plt.scatter(feat_embedded[:, 0], feat_embedded[:, 1],  feat_embedded[:, 2],c = labels_train);
	
	fig = plt.figure()
	ax = fig.add_subplot(111,projection='3d');

	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')
	plt.show();
	plt.savefig('../vis/tsne-cnn-spec-3d.png');
	plt.close()

	plt.scatter(feat_embedded[:, 2], feat_embedded[:, 0], c = labels_train);
	plt.tight_layout();
	plt.show();
	plt.savefig('../vis/tsne-cnn-spec-3d-zx.png');
	plt.close();

	plt.scatter(feat_embedded[:, 2], feat_embedded[:, 1], c = labels_train);
	plt.tight_layout();
	plt.show();
	plt.savefig('../vis/tsne-cnn-spec-3d-zy.png');
	plt.close();


# 2d
if GMM:
	feat_embedded = np.loadtxt("gmm-2d-feat_embedded.csv",dtype='float',delimiter=",");
	plt.scatter(feat_embedded[:, 0], feat_embedded[:, 1], c = labels_train);
	plt.tight_layout();
	plt.show();
	plt.savefig('../vis/tsne-gmm-spec-2d.png');

else:
	feat_embedded = np.loadtxt("gmm-3d-feat_embedded.csv",dtype='float',delimiter=",");
	plt.scatter(feat_embedded[:, 0], feat_embedded[:, 1], c = labels_train);
	plt.tight_layout();
	plt.show();
	plt.savefig('../vis/tsne-cnn-spec-2d.png');

