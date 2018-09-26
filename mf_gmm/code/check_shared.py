
directory = '/home/yichingl/mf/mf2018/mf_gmm/ctl_files';
filename1 = directory + '/' + 'sel.heliB.train_subgroups.ctl';
filename2 = directory + '/' + 'sel.heliA.test_subgroups.ctl'

with open(filename1,'r') as f1:
	with open(filename2,'r') as f2:
		for line1 in f1:
			for line2 in f2:
				if line1 == line2:
					print('found shared file: ' + line1);
