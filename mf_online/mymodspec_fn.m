
function modspecout = mymodspec(raw_audio_pylist, fs_py )

    % Converting input from python to format edible for matlab
    x = cell2mat(raw_audio_pylist);
    fs = double(fs_py);

    SR = 16000;

    x = mean(x,2);% col vector
    x = resample(x,SR,fs);
    fs = SR;

    % Resample to around 8KHz
    xlen = length(x);

    %% Construct bandpass window
    ham_t = .25; %250 ms duration window
    ham_N = floor(ham_t*fs);
    w = hamming(ham_N);
        w = w/(w'*w); % normalize energy (IMPORTANT!)
    wshift = 5; %4hz
    exp_modulator = exp(1j*2*pi*wshift/fs.*(0:ham_N-1)); %mod by 4 hz * 2pi 
    exp_modulator = exp_modulator.';
    
    w = w.*exp_modulator;

    %% Create Gammatone Filterbank
    % Make the center frequency vector
    numChannels = 128;
    lowFreq = 0; %?
    fcoefs = MakeERBFilters(fs,numChannels,lowFreq);

    % Operate on each time segment
    x_filtered =  ERBFilterBank(x, fcoefs);
    x_filtered = x_filtered.'; % flip for processing
    % time: cols, filt: rows

    % get signal envelope (envelope operates on columns)
    envt = envelope(x_filtered);
    envt = envt - mean(envt,1); % subtract average

    % bandpass signal envelope
    % filter each time vector with filter 
    bp_sig=abs(filter(w,1,envt)).^2;

    bp_sig = downsample(bp_sig,100); % downsample operates along cols

    modspecout = 10*log10(abs(bp_sig.'));

