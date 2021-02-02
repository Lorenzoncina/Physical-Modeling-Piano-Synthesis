
function output=PianoSoundSynthesis(f0,reverbType)

%PIANO HAMMER MODEL

Fs=44100; % Sampling frequency
N=65; % Number of spatial grid points
Lenght=0.62; % Length of the piano wire
MassString=3.93/1000; % Mass of the piano wire
MassHammer=2.97/1000; % Mass of the hammer
K=4.5*10^9; % Hammer stiffness coefficient
Tension=670; % Tension in the piano wire
p=2.5; % Stiffness non-linear component
alpha=0.12; % Relative striking position
b1=0.5; % Damping coefficient
b3=6.25*10^-9; % Damping coefficient
epsilon=3.82*10^-5; % String stiffness parameter
v=4; % Initial hammer velocity
R0=sqrt(Tension*MassString/Lenght); % Wave impedance of the piano wire
i0=round(alpha*N); % Striking position of the hammer
c=sqrt(Tension/(MassString/Lenght)); % Wave speed

% Defining the coefficients of the wave equation 

D=1+b1/Fs+2*b3*Fs;
r=c*N/(Fs*Lenght);
a1=( 2 - 2*r^2 + b3*Fs - 6*epsilon*N^2*r^2 )/D;
a2=( -1 + b1/Fs + 2*b3*Fs )/D;
a3=( r^2*( 1 + 4*epsilon*N^2 ) )/D;
a4=( b3*Fs - epsilon*N^2*r^2 )/D;
a5=( -b3*Fs )/D;

% Initializing some variables %

input_length=150;       %istanti temporali discreti
ys=zeros(N,input_length); % Displacement of the string
yh=zeros(1,input_length); % Displacement of the hammer
F_out=zeros(1,input_length); % Force signal output

% Initializing the values for the first few time steps of the simulation %

ys(:,1)=0;
%F(1) = 0, yh(1)=0 (forza e displ del martello 0 al primo istante temporale

%1° step (paper)
yh(2)=v/Fs;
ys(1,2)=0;   %corda fissa agli estremmi
ys(N,2)=0;
ys(2:N-1,2)=(ys(3:N,1)+ys(1:N-2,1))/2; %primo step della corda con taylor
F_out(2)=K*abs(yh(2)-ys(i0,2))^p;

%2° step (paper)
ys(1,3)=0;       %corda fissa agli estremi
ys(N,3)=0;
ys(2:N-1,3)=ys(3:N,2)+ys(1:N-2,2)-ys(2:N-1,1);
ys(i0,3)=ys(i0+1,2)+ys(i0-1,2)-ys(i0,1)+((1/Fs)^2*N*F_out(2))/MassString;
yh(3)=2*yh(2)-yh(1)-((1/Fs)^2*F_out(2))/MassHammer;
F_out(3)=K*abs(yh(3)-ys(i0,3))^p;

% Loop through the remaining time steps, implementing the finite difference
% hammer and string model.

for n=4:input_length
 ys(1,n)=0;      %estremi fissi
 ys(N,n)=0;
 ys(2,n)= a1*ys(2,n-1)+a2*ys(2,n-2)+...
 a3*(ys(3,n-1)+ys(1,n-1))+...
 a4*(ys(4,n-1)-ys(2,n-1))+...
 a5*(ys(3,n-2)+ys(1,n-2)+ys(2,n-3));
 ys(N-1,n)= a1*ys(N-1,n-1)+a2*ys(N-1,n-2)+...
 a3*(ys(N,n-1)+ys(N-2,n-1))+...
 a4*(ys(N-3,n-1)-ys(N-1,n-1))+...
 a5*(ys(N,n-2)+ys(N-2,n-2)+ys(N-1,n-3));
 ys(3:N-2,n)= a1*ys(3:N-2,n-1)+a2*ys(3:N-2,n-2)+...
 a3*(ys(4:N-1,n-1)+ys(2:N-3,n-1))+...
 a4*(ys(5:N,n-1)+ys(1:N-4,n-1))+...
 a5*(ys(4:N-1,n-2)+ys(2:N-3,n-2)+ys(3:N-2,n-3));

 ys(i0,n)= a1*ys(i0,n-1)+a2*ys(i0,n-2)+...
 a3*(ys(i0+1,n-1)+ys(i0-1,n-1))+...
 a4*(ys(i0+2,n-1)+ys(i0-2,n-1))+...
 a5*(ys(i0+1,n-2)+ys(i0-1,n-2)+ys(i0,n-3))+...
 ((1/Fs)^2*N*F_out(n-1))/MassString;

 yh(n)=2*yh(n-1)-yh(n-2)-((1/Fs)^2*F_out(n-1))/MassHammer;

 % Check for when the hammer is no longer in contact with the string 

 if (yh(n)-ys(i0,n))>0
 F_out(n)=K*abs(yh(n)-ys(i0,n))^p;
 else
 F_out(n)=0;
 end
end

% Changes the force signal into a veolcity to be fed into the digital
% waveguide model.

v=F_out/(2*R0);

%PIANO STRING MODEL

% al Loss Filter coefficient
% gl Loss Filter gain
% ad Dispersion Filter coefficient
% ap_num Number of allpass filters used in the Dispersion Filter
% offtune Variation in the Tuning Filter to make sure the three
% waveguides have different frequency
% N Length of the entire delay line of the waveguide model
% M Length of the two parallel delay lines
% P Difference between the exact delay line length reqruied and
% the actual length implemented
% C The Tuning Filter coefficient

% Initialize the output %

output_length=100000;
output=zeros(1,output_length);

% Convolves the input signal with the recorded response of the piano boday
% being knocked.

ir=audioread('AudioFile_IR\Piano_IR.wav');
v_new=conv(v,ir);
v_in=[v_new' zeros(1,length(output)-length(v_new))];



% Define/Calculate some of the parameters that will be used.

% The parameters of the filter are changed according to frequency to give a
% more consistent and normalized output.

if f0>3000
 gl=-0.997;
 ap_num=0;
 offtune=0.01;
elseif f0>1900
 gl=-0.997;
 ap_num=2;
 offtune=0.005;
elseif f0>1800
 gl=-0.997;
 ap_num=3;
 offtune=0.005;
elseif f0>1500
 gl=-0.995;
 ap_num=4;
 offtune=0.01;
elseif f0>980
 gl=-0.995;
 ap_num=6;
 offtune=0.02;
elseif f0>750
 gl=-0.993;
 ap_num=8;
 offtune=0.03;
elseif f0>390
 gl=-0.99;
 ap_num=12;
 offtune=0.04;
elseif f0>261.626
 gl=-0.985;
 ap_num=14;
 offtune=0.06;
elseif f0>200
 gl=-0.98;
 ap_num=16;
 offtune=0.09;
elseif f0>150
 gl=-0.975;
 ap_num=18;
 offtune=0.13;
elseif f0>120
 gl=-0.968;
 ap_num=20;
 offtune=0.18;
else
 gl=-0.96;
 ap_num=20;
 offtune=0.25;
end
al=-0.001;
ad=-0.30;
%lunghrzza esatta della delay line
N_exact=((2*pi+ap_num*atan(((ad^2-1)*sin(2*pi*f0/Fs))/...
     (2*ad+(ad^2+1)*cos(2*pi*f0/Fs))))/(2*pi*f0/Fs));
M=floor(N_exact/2);
P=N_exact-2*M;
C=(1-P)/(1+P);
i0=round(alpha*M);


% Defines the transfer function for the delays and filters used:
%
% DL1 The delay line representing the segment of the string from the
% agraffe to the point of contact with the hammer
% DL2 The delay line representing the segment of the string from the
% point of contact with the hammer to the bridge
% Hl The Loss Filter
% Hd One of the the allpass filters that make up the Dipersion Filter
% Hfd Hfd1, Hfd2 and Hfd2 are the Tuning Filters used to tune the
% fundamental frequency. The 3 filters are each used in one of the
% 3 parallelly connected digital waveguide models.
%
% H H1, H2 and H3 are the Reflections Filters for each of the 3
% digital waveguide models. They are made up of the Loss Filter,
% the Dipersion Filter and the Tuning Filter.
%

z=tf('z',1/Fs);
DL1=(z)^-(M-i0);
DL2=(z)^-(i0);
Hl=gl*(1+al)/(1+al*z^-1);
Hd=(ad+z^-1)/(1+ad*z^-1);
%le tre waveguide in paralelo sono uguali in tutto e per tutto,
%ma hanno 3 differenti tuning filter, in quanto modellano il piano reale
Hfd1=(C+z^-1)/(1+C*z^-1);
Hfd2=(C*(1+offtune)+z^-1)/(1+C*(1+offtune)*z^-1);
Hfd3=(C*(1-offtune)+z^-1)/(1+C*(1-offtune)*z^-1);
if C*(1+offtune)>=1 % Makes sure the coefficient of the
 Hfd2=(1+z^-1)/(1+1*z^-1); % Tuning Filter does not excede 1
end
%costituisco i 3 reflection filter (il filtro totale è dato dal prodtto
%delle 3 funzioni di trasferimento di ciascuno (per comporre LTI in serie
%basta moltiplicare le funzioni di trasferimento
H1=Hl*Hd^ap_num*Hfd1;
H2=Hl*Hd^ap_num*Hfd2;
H3=Hl*Hd^ap_num*Hfd3;




% The filters are then combined according to the digital waveguide model of
% the piano string.

DW1=DL1/(1+H1*DL1*DL1*DL2*DL2)+DL2*DL2*DL1*(-1)/(1+H1*DL1*DL1*DL2*DL2);
DW2=DL1/(1+H2*DL1*DL1*DL2*DL2)+DL2*DL2*DL1*(-1)/(1+H2*DL1*DL1*DL2*DL2);
DW3=DL1/(1+H3*DL1*DL1*DL2*DL2)+DL2*DL2*DL1*(-1)/(1+H3*DL1*DL1*DL2*DL2);

% The Digital Waveguide filter is then used on the input velocity.

[b,a]=tfdata(DW1,'v');
output1=filter(b,a,v_in);
[b,a]=tfdata(DW2,'v');
output2=filter(b,a,v_in);
[b,a]=tfdata(DW3,'v');
output3=filter(b,a,v_in);


% Output of the three digital waveguides re summed together. The sum is
% then normalized. The final output is then played.

output=output1+output2+output3;
output=output/max(abs(output))*(1 - 1/32768);
%soundsc(output,Fs)



%%%REVERB CHOICE
switch reverbType
    case 'No riverbero'
        soundsc(output,Fs);
        %audiowrite('./output.wav', output, Fs)
    case 'Riverbero 1'
        %aula conferenze (lecture)
        airpar.rir_type = 1;
        airpar.room = 4;
        airpar.fs = Fs;
        airpar.channel = 1;
        airpar.head = 1;
        airpar.rir_no = 6;
       
        [h_air,air_info] = load_air(airpar);
        
        %implemento convoluzione (implementare funzione da richiare) 
        %funzione riverberazione
        %reverbOutput = conv(output, h_air);
        reverbOutput = reverb(output, h_air,Fs, f0);
        soundsc(reverbOutput,Fs)
    case 'Riverbero 2'
        %stairway
        airpar.rir_type = 1;
        airpar.fs = Fs;
        airpar.room = 5;
        airpar.channel = 1;
        airpar.rir_no = 2;
        airpar.azimuth = 15;
        airpar.head = 1;
        [h_air,air_info] = load_air(airpar);
        
        reverbOutput = reverb(output, h_air,Fs, f0);
        soundsc(reverbOutput,Fs)
    case 'Riverbero 3'
        airpar.rir_type = 1;
        airpar.fs = Fs;
        airpar.room = 2; %stanza ufficio
        airpar.channel = 1;
        airpar.rir_no = 2;
        airpar.head = 1;
        [h_air,air_info] = load_air(airpar);
        
        reverbOutput = reverb(output, h_air,Fs, f0);
        soundsc(reverbOutput,Fs)

end

%visualizza o meno il plot della stringa
plot = false;

if plot
    stringPlot(M,output_length, [-.1 -.25 -.1] , alpha );
end


%%%%%%%%%% reverb
% airpar.rir_type = 1;
% airpar.room = 4;
% airpar.fs = Fs;
% airpar.channel = 1;
% airpar.head = 1;
% airpar.rir_no = 5;
% [h_air,air_info] = load_air(airpar);
% 
% reverbOutput = conv(output, h_air);
% soundsc(reverbOutput,Fs)

% %Plot the string's movement
% if(enablePlot)
% %     ii = 0:(M-1);
% %     for t = 1:output_length
% %         plot(ii, output);
% %     end
%     figure();
%     subplot(2,1,1);
%     plot(output);
%       
% end




