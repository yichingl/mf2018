% matlab function that gets called by testpy and testpy2.py

function a = test(v)

	%for inputting scalar
	a = double(v) ~= 5

	% for inputting an np.array (convert to list first, pass in as cell)
	%a = cell2mat(vector);
	%a = sum(a);

	% for testpy
	%a = [vector; zeros(1,length(vector))];