function y = stringPlot(len,d,r,x)

%    len is the length of each waveguide in samples.  
%    d number of output samples to produce (default 10000).
%    r  refelection filter 
%    x is optional initial shape of the string (triangular di default)
%    plot if present and nonzero plots waves & string at each sample time.


% where to read output from - residual at bridge in this case
pickup = len;

firlen = length(r);
firHlen = floor(firlen/2);


% Make right-going rail long enough to hold [r] points centered 
% at the end point (for calculating zero-phase reflection)
rlen = len + firlen - firHlen;
right = zeros(1,rlen);

left = zeros(1,len);

% Does x specify pluck point?
if length(x) == 1
  % Interpret as proportion along string at which it is plucked
  pluck = x*(len-1);
  x = [ [0:floor(pluck)]/pluck, ...
	(len - 1 - [(floor(pluck)+1):(len-1)])/(len - 1 - pluck) ];
end

% Initialization
if length(x) < len
  dl = len - length(x);
  x = [zeros(1,floor(dl/2)), x, zeros(1,ceil(dl/2))];
end

% Because initial velocity profile is flat, initial displacement 
% profile is equal in leftgoing and rightgoing waves.
left(1:len) = x(1:len)/2;
right(1:len) = x(1:len)/2;
% fill in extra point for bridge filter
right(len+1) = 0;

% Initialize output
y = zeros(1,d);

% Initialize variables for display
pkval = max(abs(x));
ii = 0:(len-1);

% Execute waveguide
for t = 1:d
  

   % Plot left and right-moving waves, and their sum
   plot(ii, left, ii, right(1+ii), ii, left+right(1+ii));
   % Make sure the axis scaling stays the same for each plot
   axis([0 len-1 -pkval pkval]);
   pause(0.000005);
 
    
  % Move left-hand moving wave one step left; append dummy value for now
  left = [left(2:len),0];
  % At 'nut' (left-hand end), assume perfect reflection, so new value 
  % of right-moving wave is negative of new value at nut of left-moving
  nut = -left(1);
  % Move right-moving wave one step (including extra point off end)
  right = [nut, right(1:len)];
  % Apply 'bridge' filter to end points of right-moving wave to get 
  % new value to fill in to end of left-moving wave
  % One point of convolution
  bridge = r * right( (len-firHlen-1)+[1:firlen] )';
  
  %bridge = filter(r,den, right( (len-firHlen-1)+[1:firlen] ));
  disp(bridge);
  left(len) = bridge;  

  % Read output
  y(t) = left(pickup) + right(pickup);

end