import re; #for extracting substrings

# code for creating group names for TRAINING data

# USER INPUT
group_name = 'allcat_bal';
type_name = 'heli'; #'heli'
file_name = 'heli'; #boat #heli #boatA #heliA
test_or_train = 'train';

boat_starts = ['boat/'+test_or_train+'/','extra/boat/'];
heli_starts = ['heli/'+test_or_train+'/','extra/heli/'];
end = '_'; # underscore between subgroup name and number



ctl_path = '../ctl_files/' + group_name +'.' + file_name+'.'+test_or_train+'.ctl';
subgroups_path = '../ctl_files/' + group_name +'.' + file_name+'.'+test_or_train+'_subgroups.ctl';

# clear file
with open(subgroups_path,'w') as s:
	print ;

with open(ctl_path) as f: 
	with open(subgroups_path,'a') as s:

		# extract file name
		for ctl_name in f:
			print(ctl_name);
			# process boat files, extract from before underscore
			if type_name == 'boat':
				for boat_start in boat_starts:
					results1 = re.search(boat_start+'(.*)'+end, ctl_name);
					if results1 != None:
						subgroup_name = results1.group(1);

			# process heli files, extract from before underscore
			if type_name == 'heli':
				for heli_start in heli_starts:
					results1 = re.search(heli_start+'(.*)'+end, ctl_name);
					if results1 != None:
						subgroup_name = results1.group(1);
			# look for second underscore
			results2 = re.search('(.*)'+end, subgroup_name);

			# if found second underscore, use new result
			if results2 != None:
				subgroup_name = results2.group(1);
			s.write(subgroup_name+'\n');
