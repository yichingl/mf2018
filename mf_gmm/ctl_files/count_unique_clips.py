from collections import Counter
import numpy as np


filenames = ['allcat_bal.boat.train_subgroups.ctl',
	'allcat_bal.heli.train_subgroups.ctl', 
	'allcat_bal.boat.test_subgroups.ctl', 'allcat_bal.heli.test_subgroups.ctl'];

for filename in filenames:
	with open(filename) as f:
	    clip_names = f.read().splitlines()

	print('===============')
	print('Info for: '+filename);
	print(Counter(clip_names).keys()) # equals to list(set(words))
	print(Counter(clip_names).values()) # counts the elements' frequency
	print('===============\n')