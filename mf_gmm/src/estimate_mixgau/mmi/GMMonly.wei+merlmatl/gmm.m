function [g,p] = gmm( x, k, a, m, r, st, pl, thr)
% Gaussian mixture model
%  function [g,p] = gmm( x, k, a, m, r, st, pl, thr)
%
% Training:
%  x is the input data (dimensions x samples)
%  k (or k.k) is the number of Gaussians to learn
%  a (or k.a) are the initial Gaussian priors (default = ones(1, k)/k)
%  m (or k.m) are the initial Gaussian means (default from kmeans)
%  r (or k.r) are the initial Gaussian covariances (default = eye(k))
%  st (or k.st) are the covariance forms ('f','d','u'=full/diagonal/unit, default = 'f')
%  pl (or k.pl) is 1 to turn on plotting info, or 0 otherwise (default = 0)
%  thr (or k.thr) is the log likelihood threshold to stop training (default = .001)
%  g is the Gaussian mixture data structure
%
% Classification:
%  x is the input data
%  k is the Gaussian mixture data structure
%  g are the Gaussian assignments to the samples
%  p is the log probability of the model
%
%  Paris Smaragdis, 14-May-03

% Get traininng parameters if we have any in a struct
p = 0;
if (nargin == 2) & isstruct( k) & ~isfield( k, 'fun')
	if isfield( k, 'a')
		a = k.a;
	else
		a = [];
	end
	if isfield( k, 'm')
		m = k.m;
	else
		m = [];
	end
	if isfield( k, 'r')
		r = k.r;
	else
		r = [];
	end
	if isfield( k, 'st')
		st = k.st;
	else
		st = 'f';
	end
	if isfield( k, 'thr')
		thr = k.thr;
	else
		thr = .001;
	end
	if isfield( k, 'pl')
		pl = k.pl;
	else
		pl = 0;
	end
	k = k.k;
end

% Training
if ~isstruct( k)
	% Do default values
	if nargin < 3
		a = [];
	end
	if nargin < 4
		m = [];
	end
	if nargin < 5
		r = [];
	end
	if nargin < 6
		st = 'f';
	end
	if nargin < 7
		pl = 0;
	end
	if nargin < 8
		thr = .001;
	end

	% and train
	g = gmm_train( x, k, a, m, r, st, pl, thr);

% Evaluation
else
	[g,p] = gmm_run( x, k);

end


%-------------------------------
% Gaussian mixtures training
function g = gmm_train( x, k, a, m, r, st, plt, thr)

% Setup some plotting
if plt
	set( gcf, 'DoubleBuff', 'on');
end

[d,n] = size( x);

% If we don't provide weights, assume equality
if isempty( a)
	a = ones( k, 1)/k;
end

% Likewise for means, use k-means to figure an initial values
if isempty( m)
if 0
	[km,cl] = pkmeans( x, k, 100);
	for i = 1:k
		m{i} = km.ctr(:,i);
	end
else
	for i = 1:k
		m{i} = x(:,i);
	end	
end
end

% And likewise for covariances, start with diagonals
if isempty( r)
	globalcov = inv(pov(x));
	lgc = log(det(globalcov));
	for i = 1:k
		ir{i} = globalcov;
		logdir(i) = lgc;
	end
end
ll = zeros( k, n);
lk = NaN;

% Start iterating
for ep = 1:200:-1
	% Expectation
	for i = 1:k
		xm = repop( x, m{i}, '-');
		switch st
			case {'f' , 'sf'}
				dist = sum( xm .* (ir{i}*xm));
			case {'d', 'sd'}
				dist = sum( xm .* repop( xm, diag( ir{i}), '*'));
			case 'u'
				dist = sum( xm .* xm);
		end
		ll(i,:) = -.5*d*log(2*pi) + .5*logdir(i) - dist/2;
	end
	ll = repop( ll, log(a), '+');
	mll = max(ll);
	ll = ll - repmat( mll, k, 1);
	fl = lsum( ll')';
	pl = ll - repmat( fl, k, 1);
	lk(end+1) = sum( fl(:)) + sum(mll);
	p = exp( pl);

	% Maximization
	ss = 0;
	for i = 1:k
		ps = sum( p(i,:));

		% Weights
		a(i) = ps/n;

		% Means
		m{i} = x * p(i,:)' / ps;

		% Variances
		xm = repop( x, m{i}, '-');
		pp = repmat( p(i,:), d, 1);
		switch st
			case 'f'
				ss = (pp.*xm) * xm' / ps;
				ir{i} = inv( ss);
				logdir(i) = log( det( ir{i}));
			case 'd'
				ss = sum( pp.*xm.^2, 2) / ps;
				ir{i} = diag( 1./ss);
				logdir(i) = sum( log( diag( ir{i})));
			case 'sf'
				ss = ss + (pp.*xm) * xm';
			case 'sd'
				ss = ss + sum( pp.*xm.^2, 2);
		end
	end

	% If covariances are shared
	if st(1) == 's'
		ss = ss/n;
		switch st
			case 'sf'
				ss = inv(ss);
				logdetss = log(det(ss));
			case 'sd'
				ss = diag(1./ss);
				logdetss = sum(log(diag(ss)));
		end
		for i = 1:k
			ir{i} = ss;
			logdir(i) = logdetss;
		end
	end

	% Some plotting
	if plt == 1
		if d < 4
			[ll,cl] = max( p);
			subplot( 2, 1, 1), scatplot( x, cl);
			subplot( 2, 1, 2)
		end
		if ep > 10
			plot( ep-10:ep, lk(end-10:end)); grid on
		else
			plot( lk); grid on
		end
		title( sprintf( 'Log Likelihood, Last delta: %f', lk(end)-lk(end-1)'));
		drawnow
    end
    disp(['ep = ' num2str(ep) ' : lkd = ' num2str(lk(end))]);
	% Time to quit?
	if abs( lk(end)-lk(end-1)) < thr
		break
	end
end

% Lemme know if we never made it
if ep == 200
	warning( 'Log likelihood threshhold not achieved after 200 iterations');
end

% Save to a data structure
g.a = a;
g.m = m;
for i = 1:k
	g.ir{i} = ir{i};
	g.gc(i) = -.5*rows(x)*log(2*pi) + .5*log( det( g.ir{i}));
end
g.fun = 'gmm';
g.t = st;


%-------------------------------
% Gaussian mixtures evaluation
function [c,p] = gmm_run( x, g)

p = zeros( length( g.a), cols( x));
switch g.t
	% Full covariances
	case {'f' 'sf'}
		for i = 1:length( g.a)
			xm = repop( x, g.m{i}, '-');
			d = sum( xm .* (g.ir{i}*xm), 1);
			p(i,:) = log( g.a(i)) - d/2 + g.gc(i);
		end
	% Diagonal covariances
	case {'d', 'sd'}
		for i = 1:length( g.a)
			xm = repop( x, g.m{i}, '-');
			d = sum( repop( xm.*xm, diag( g.ir{i}), '*'), 1);
			p(i,:) = log( g.a(i)) - d/2 + g.gc(i);
		end
	% Unit covariances
	case 'u'
		for i = 1:length( g.a)
			xm = repop( x, g.m{i}, '-');
			d = sum( xm .* xm, 1);
			p(i,:) = log( g.a(i)) - d/2 + g.gc(i);
		end
end
[c,c] = max( p,[],1);
p = lsum( p')';
