%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                           %
%      LEAST MEAN SQUARES (LMS) ALGORITHM & EQUALIZATION    %
%      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   %
%                                                           %
% AUTHOR         :   SRIDHAR KRISHNAN                       %
% STUDENT ID # :     12345                                  %
% File           :   equ.m                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
%% Assumptions for Channel Capacity W, No. of Samples, No. of trials %%
W = 2.9;                    % Channel Capacity
alp = 2;                    % alpha
forder = 7;                 % filter order
tsamples = 700;             % total # of samples
delay = 5;                  % delay for desired response
trials = 10;                % No. of independent trials
de = zeros(1,tsamples);     % Intializing the squared error
nvari = 0.001;              % Noise variance
dvari = 1;                  % Input data variance
taps = zeros(1,forder);     % Initializing the tap weight counter
muu = 0.0016;                 % Fill the appropriate mu
errorsq = 0;                % Mean squared Error
op = zeros(1,forder);
e = zeros(1,forder);

for trys = 1:trials           % Loop for No. of trials
    u = zeros(1,forder-1);     % intialization of tap I/p for time < 1
   a = zeros(1,forder-1);     % intialization of data i/p for time < 1
    for j = 1:forder
        w(forder,j) = 0;      % intializing tapweight vector w(0) = 0
    end
%%%%%%%%%% MODULE 1 %%%%%%%%%%%%%%%%
for sampgen = forder : tsamples+forder   %% Generates Random data & desired resp
   a(sampgen) = fix(rand+0.5)*2-1;
end
 
    switch trys
        case 1
            W = 2.9;
        case 2
            W = 3.3;
        case 3
            W = 3.6;
        case 4
            delay = 9;
        case 5
            nvari = 0.001;
        case 6
            nvari = 0.004;
        case 7
            muu = 0.2;
        case 8
            muu = 0.16;
        case 9
            muu = 0.11;
        case 10
            muu = 0.09;
    end
            

d = [zeros(1,delay) a];
for sampgen = forder : tsamples+forder %% Generates Random noise
 v(sampgen) = fix(rand+0.5)*2*sqrt(nvari)-sqrt(nvari);
end
%% Impulse Response of the Channel for given W %%
                                                 
for k = 1:3
h(k) = 1/2*(1+cos(2*pi/W*(k-2)));
end;
%% Mu
%%muu = 0.15;   %%% Fill the appropriate mu

%%%%%%% MODULE 2 : UPDATING THE WEIGHTS & PLOTTING THE LEARNING CURVE %%%%%%%
%% Output from Channel
y = conv2(a,h,'same');
%% Input to the equalizer (Channel o/p + Random Noise)
eq = y + v;

%% desired response
adelay = a(delay:end);

for n = forder : size(adelay,2)
    u = eq(n:-1:n-forder+1) ;
    %% Finding the Estimate
    op(n)= taps * u';
    
%% Calculating the Error (desired resp - estimate)
    e(n) = adelay(n) - op(n) ;
    
%% Updating the tap weights
    taps = taps + muu * u * e(n) ;
end 
%%Checking:

%% Equalizer Output with updated tap weights
eqFinal = conv2(eq,taps,'same');
%% Combined Response of Channel and Equalizer
hCombined = conv(h,taps);
%% Calculating the error square for each trial & adding it up
errorsq  = errorsq + e.^2;
end %%% Loop for number of trials (try variable)
%% Calculating the Mean square error (ie total err at instant i / # of trials)
errorsq = errorsq/10;
%% Plotting the Learning Curve
semilogy(errorsq);
title('Learning Curve') ;
xlabel('Samples')
ylabel('Mean Squared Error')

%%%%%   Plotting  the channel responses %%%%
plot(h);
title(sprintf('Channel response for W = %d', W)) ;
xlabel('Samples')
ylabel('Channel Response')
%%%%% Plotting the equalizer responses   %%%%
plot(taps);
title('Equalizer Response') ;
xlabel('Samples')
ylabel('Equalizer Response')
%%%%% Plotting   the combined responses %%%%%

plot(hCombined);
title('Combined Response') ;
xlabel('Samples')
ylabel('Combined (Equalizer and Channel) Response')