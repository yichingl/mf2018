import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

labels_train = np.loadtxt("labels.csv",dtype='float');

GMM=0;
# 3d
if GMM:
	fig = plt.figure()
	ax = fig.add_subplot(131,projection='3d');
	
	feat_embedded = np.loadtxt("gmm-3d-feat_embedded.csv",dtype='float',delimiter=",");

	xs = np.asarray(feat_embedded[:,0]);
	ys = np.asarray(feat_embedded[:,1]);
	zs = np.asarray(feat_embedded[:,2]);
	ax.scatter(xs,ys,zs,c=labels_train);
	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')

	rad_split = 13;

	#fig = plt.figure()
	ax = fig.add_subplot(132,projection='3d');
	radius = np.sqrt(np.square(xs)+np.square(ys)+np.square(zs));
	print(radius[30:40]);
	rad_ind = radius<=rad_split;
	selected_feat = feat_embedded[rad_ind,:];
	print(rad_ind[30:40]);
	selected_labels = labels_train[rad_ind];
	print(selected_feat[0:10])
	xs_sel = selected_feat[:,0];
	ys_sel = selected_feat[:,1];
	zs_sel = selected_feat[:,2];
	ax.scatter(xs_sel,ys_sel,zs_sel,c=selected_labels);
	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')


	#fig = plt.figure()
	ax = fig.add_subplot(133,projection='3d');
	radius = np.sqrt(np.square(xs)+np.square(ys)+np.square(zs));
	print(radius[30:40]);
	rad_ind = radius>rad_split;
	selected_feat = feat_embedded[rad_ind,:];
	print(rad_ind[30:40]);
	selected_labels = labels_train[rad_ind];
	print(selected_feat[0:10])
	xs_sel = selected_feat[:,0];
	ys_sel = selected_feat[:,1];
	zs_sel = selected_feat[:,2];
	ax.scatter(xs_sel,ys_sel,zs_sel,c=selected_labels);
	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')
	plt.show();
	plt.close()
else:	
	
	fig = plt.figure()
	ax = fig.add_subplot(131,projection='3d');
	
	feat_embedded = np.loadtxt("cnn-3d-feat_embedded.csv",dtype='float',delimiter=",");

	xs = np.asarray(feat_embedded[:,0]);
	ys = np.asarray(feat_embedded[:,1]);
	zs = np.asarray(feat_embedded[:,2]);
	ax.scatter(xs,ys,zs,c=labels_train);
	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')

	rad_split = 13;

	#fig = plt.figure()
	ax = fig.add_subplot(132,projection='3d');
	radius = np.sqrt(np.square(xs)+np.square(ys)+np.square(zs));
	print(radius[30:40]);
	rad_ind = radius<=rad_split;
	selected_feat = feat_embedded[rad_ind,:];
	print(rad_ind[30:40]);
	selected_labels = labels_train[rad_ind];
	print(selected_feat[0:10])
	xs_sel = selected_feat[:,0];
	ys_sel = selected_feat[:,1];
	zs_sel = selected_feat[:,2];
	ax.scatter(xs_sel,ys_sel,zs_sel,c=selected_labels);
	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')


	#fig = plt.figure()
	ax = fig.add_subplot(133,projection='3d');
	radius = np.sqrt(np.square(xs)+np.square(ys)+np.square(zs));
	print(radius[30:40]);
	rad_ind = radius>rad_split;
	selected_feat = feat_embedded[rad_ind,:];
	print(rad_ind[30:40]);
	selected_labels = labels_train[rad_ind];
	print(selected_feat[0:10])
	xs_sel = selected_feat[:,0];
	ys_sel = selected_feat[:,1];
	zs_sel = selected_feat[:,2];
	ax.scatter(xs_sel,ys_sel,zs_sel,c=selected_labels);
	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')
	plt.show();
	plt.close()




	#plt.scatter(feat_embedded[:, 2], feat_embedded[:, 0], c = labels_train);
	#plt.tight_layout();
	#plt.show();
	#plt.savefig('../vis/tsne-cnn-spec-3d-zx.png');
	#plt.close();

	#plt.scatter(feat_embedded[:, 2], feat_embedded[:, 1], c = labels_train);
	#plt.tight_layout();
	#plt.show();
	#plt.savefig('../vis/tsne-cnn-spec-3d-zy.png');
	#plt.close();


# 2d
if GMM:
	feat_embedded = np.loadtxt("gmm-2d-feat_embedded.csv",dtype='float',delimiter=",");
	plt.scatter(feat_embedded[:, 0], feat_embedded[:, 1], c = labels_train);
	plt.tight_layout();
	#plt.show();
	plt.savefig('../vis/tsne-gmm-spec-2d.png');

else:
	feat_embedded = np.loadtxt("gmm-3d-feat_embedded.csv",dtype='float',delimiter=",");
	plt.scatter(feat_embedded[:, 0], feat_embedded[:, 1], c = labels_train);
	plt.tight_layout();
	#plt.show();
	plt.savefig('../vis/tsne-cnn-spec-2d.png');

