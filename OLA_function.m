function outputArg1 = OLA(input_signal,impulse_response, Fs, f0)

% L parameter (window-related)
L = 4;

% length of the analysis window
M = ceil(L * Fs / f0);

% analysis window
win = hamming(M);

% analysis hop size
R = round(M / 4);

%fare uno switch con frequenze sopra o sotto 600 Hz(vedi slide sarti
%sinusoidal analisys)

%JND just noticeable difference
JND = 3;

% minimum fft length to be under JND (better explained in report)
min_fft_length = Fs / (2 * JND);

N_fft = 2^(ceil(log2(min_fft_length)));
%N_fft_h=2^(ceil(lenght(impulse_response)))

%number of frames
n_frames = floor(((length(input_signal) - N_fft)/R))+1;


%filter transform
H=fft(impulse_response,N_fft);

%Nyquist
N_nyq = N_fft/2;
% 

%length of block
Lconv = M + length(impulse_response) - 1;

%total length of output
y_oa = zeros(1, Lconv*num_frames);

for m = 0:n_frames-1
% signal windowing 
    x_m = x(m*R + 1:m*R + M);
    y_m = win.*x_m;
        
    % compute the FFT
    X = fft(y_m, N_fft);
    
    % retain only the positive frequencies up to the Nyquist index
    X = X(1:N_nyq);
    
    
   
    y_block = ifft(X .* H);
    
    y_oa(1 + m *(M):m * M + Lconv) = y_oa(1 + m *(M): m * M + Lconv) + y_block;  
    
end

outputArg1=y_oa;

end
