import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

labels_train = np.loadtxt("labels.csv",dtype='float');

# Group name
group_name = 'allcat_bal';
feat_type = 'spec';
n_clusters = 16;

if feat_type == 'spec':
	n_feats = 133;
	n_timeframes = 669;

elif feat_type == 'mel':
	n_feats = 128;
	n_timeframes = 173;

elif feat_type == 'mod':
	n_feats = 128;
	n_timeframes = 320;
full_name = group_name+'.'+feat_type+'.'+str(n_clusters)+'n_cl'+'.'+str(n_feats)+'n_feats';
name = full_name+'.TSNE';

GMM=1;
# 3d
if GMM:
	fig = plt.figure()
	ax = fig.add_subplot(131,projection='3d');
	
	feat_embedded = np.loadtxt(name+".3D.csv",dtype='float',delimiter=",");

	xs = np.asarray(feat_embedded[:,0]);
	ys = np.asarray(feat_embedded[:,1]);
	zs = np.asarray(feat_embedded[:,2]);
	ax.scatter(xs,ys,zs,c=labels_train);
	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')

	rad_split = 10;

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
	feat_embedded = np.loadtxt(name+".2D.csv",dtype='float',delimiter=",");
	plt.scatter(feat_embedded[:, 0], feat_embedded[:, 1], c = labels_train);
	plt.tight_layout();
	plt.savefig(name+".2D.jpg");
	plt.show();

else:
	feat_embedded = np.loadtxt("gmm-3d-feat_embedded.csv",dtype='float',delimiter=",");
	plt.scatter(feat_embedded[:, 0], feat_embedded[:, 1], c = labels_train);
	plt.tight_layout();
	plt.savefig('../vis/tsne-cnn-spec-2d.png');
	plt.show();
	

