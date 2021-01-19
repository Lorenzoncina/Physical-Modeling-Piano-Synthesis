function outputArg1 = reverb(input_signal,impulse_response, Fs, f0)

% shape factor parameter (window-related)
 shapeFact= 2;

% length of the analysis window (block)
L = ceil(shapeFact * Fs / f0);

Lconv = L + length(impulse_response) - 1;

% how many non-overlapping blocks of length L are there in input_signal? 
num_blocks = ceil(length(input_signal) / L);

if length(input_signal) < num_blocks * L
    % pad x with zeros 
    x_pad = padarray(input_signal, [0, num_blocks*L - length(input_signal)], 'post');
else
    x_pad = input_signal;
end

y_oa = zeros(1, Lconv*num_blocks);

reverb = fft([impulse_response, zeros(1, Lconv- length(impulse_response))]);

for b = 1:num_blocks
    
    block_f = fft([x_pad(1 + (b-1)*(L):(b-1)*L + L), zeros(1, Lconv-L)]);
    
    y_block = ifft(block_f .* reverb);
    y_oa(1 + (b-1)*(L):(b-1)*L + Lconv) = y_oa(1 + (b-1)*(L):(b-1)*L + Lconv) + y_block;   
    
end

y_oa = y_oa(1:length(input_signal) + length(impulse_response) -1);
outputArg1=y_oa;

end
