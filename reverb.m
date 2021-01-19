function outputArg1 = reverb(input_signal,impulse_response, Fs, f0)

% shape factor parameter (window-related)
L = 4;

% length of the analysis window
M = ceil(L * Fs / f0);

% analysis window
win = hamming(M);

% analysis hop size
R = round((M-1)/2);

%JND just noticeable difference
if(f0<600)
    JND = 3;
else
    JND = (f0 / 100)*0.5;
end

% minimum fft length to be under JND 
min_fft_length = Fs / (2 * JND);
N_fft = 2^(ceil(log2(min_fft_length)));

%number of frames
n_frames = floor(((length(input_signal) - N_fft)/R))+1;

%length of a convolved block
Lconv = M + length(impulse_response) - 1;

%zero padding ir
%number_of_zeros = N_fft-M;
% ir = [impulse_response zeros(1,number_of_zeros)];
%padarray(impulse_response, [0, number_of_zeros], 'post');
%transfer function
H=fft(impulse_response,N_fft);

%Nyquist frequency
N_nyq = N_fft/2;

%reverbered signal
y_oa = zeros(1, Lconv * n_frames);

for m = 0:n_frames-1
% signal windowing 
    x_m = input_signal(m*R + 1:m*R + M);
    y_m = win.*x_m;
    
    %zero padding
    %y_pad = [y_m zeros(1,number_of_zeros)];
    %padarray(y_m, [0, number_of_zeros], 'post');
    % compute the FFT
    X = fft(y_m, N_fft);
    %plot(X);
    % retain only the positive frequencies up to the Nyquist index
    %X = X(1:N_nyq);
    disp(size(X));
    disp(size(H));
    y_block = ifft(X' .* H);
    
    %overlapp and add method
    y_oa(1 + m *(M):m * M + Lconv) = y_oa(1 + m *(M): m * M + Lconv) + y_block;  
    
end

outputArg1=y_oa;

end
