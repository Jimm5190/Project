% MATLAB script for Illustrative Problem 7.10.
clear all;
close all;

SNRindB1=0:2:15;
SNRindB2=0:0.1:15;
SNRindB3=0:2:15;
M=16;
k=log2(M);
for i=1:length(SNRindB1),
    smld_err_prb_3d(i)=cm_sm41(SNRindB1(i));	% simulated error rate
end;

for i=1:length(SNRindB3),
    smld_err_prb_4d(i)=cm_sm42(SNRindB3(i));	% simulated error rate
end;

for i=1:length(SNRindB2),
    SNR=exp(SNRindB2(i)*log(10)/10);    	% signal-to-noise ratio
    % theoretical symbol error rate
    theo_err_prb(i)=4*qfunc(sqrt(3*k*SNR/(M-1)));
    %theo_err_prb(i)=4*erfc(sqrt(3*k*SNR/(M-1)/2))/2;
end;

% Plotting commands follow.
semilogy(SNRindB1, smld_err_prb_3d, '*', 'color', 'b');
hold on;
semilogy(SNRindB3, smld_err_prb_4d, '*', 'color', 'r');
semilogy(SNRindB2,theo_err_prb,'Color','g');
hold off;
legend('Simulated symbol error rate (3A)', 'Simulated symbol error rate (4A)', 'Theoretical symbol error rate');
xlabel('E_b/N_0 in dB','fontsize',16,'fontname','Helvetica');
ylabel('Error Probability','fontsize',16,'fontname','Helvetica');
title('Performance of a 16-QAM system from Monte Carlo simulation','fontsize',12,'fontname','Helvetica');
fname = 'prob_07_10(411086040)_add_4d.png';
print (fname, '-dpng'); 

function [p]=cm_sm41(snr_in_dB)
% [p]=cm_sm41(snr_in_dB)
%		CM_SM41  finds the probability of error for the given
%   		value of snr_in_dB, SNR in dB.
N=10000;
d=1;				  	% min. distance between symbols
Eav=10*d^2;		 	  	% energy per symbol
snr=10^(snr_in_dB/10);	 	  	% SNR per bit (given)
sgma=sqrt(Eav/(8*snr));	  	  	% noise variance
M=16;
% Generation of the data source follows.
for i=1:N
  temp=rand;		        	  	% a uniform R.V. between 0 and 1
  dsource(i)=1+floor(M*temp);	  	% a number between 1 and 16, uniform 
end
% Mapping to the signal constellation follows.
mapping=[-3*d 3*d;
	      -d  3*d;
           d  3*d;
	     3*d  3*d;
	      -3*d  d;
	        -d  d;
             d  d;
	       3*d  d;
     	 -3*d  -d; 
 	       -d  -d; 
	        d  -d;
          3*d  -d;
   	   -3*d  -3*d;
	     -d  -3*d;
	      d  -3*d;
	    3*d  -3*d];
for i=1:N
  qam_sig(i,:)=mapping(dsource(i),:);
end
% received signal
for i=1:N
  n  = [sgma*randn sgma*randn];
  %[n(1) n(2)]=gngauss(sgma);
  r(i,:) = qam_sig(i,:) + n;
end
% detection and error probability calculation
numoferr=0;
for i=1:N
    
    % write your decoder here
    % The decision is named decis
    % Hint.
    % Given received signal vector r(i,:),
    % find the minimum distance among r(i,:) and the 16 signal points in
    % mapping(j,:), j=1,2,... 16. If the distance to mapping(k,:) is
    % minimum, decis=k
    min_distance = inf; % Initialize minimum distance to a large value
    for k = 1:16
        % Calculate Euclidean distance between received signal and each point in mapping
        distance = norm(r(i,:) - mapping(k,:));
        if distance < min_distance
            decis = k; % Update decision if current distance is smaller than previous minimum
            min_distance = distance;
        end
    end
  if (decis~=dsource(i))
    numoferr=numoferr+1;
  end
end
p=numoferr/(N);
end

function [p]=cm_sm42(snr_in_dB)
% [p]=cm_sm41(snr_in_dB)
%		CM_SM41  finds the probability of error for the given
%   		value of snr_in_dB, SNR in dB.
N=10000;
d=1;				  	% min. distance between symbols
Eav=10*d^2;		 	  	% energy per symbol
snr=10^(snr_in_dB/10);	 	  	% SNR per bit (given)
sgma=sqrt(Eav/(8*snr));	  	  	% noise variance
M=16;
% Generation of the data source follows.
for i=1:N
  temp=rand;		        	  	% a uniform R.V. between 0 and 1
  dsource(i)=1+floor(M*temp);	  	% a number between 1 and 16, uniform 
end
% Mapping to the signal constellation follows.
mapping=[-4*d 4*d;
	      -d  4*d;
           d  4*d;
	     4*d  4*d;
	      -4*d  d;
	        -d  d;
             d  d;
	       4*d  d;
     	 -4*d  -d; 
 	       -d  -d; 
	        d  -d;
          4*d  -d;
   	   -4*d  -4*d;
	     -d  -4*d;
	      d  -4*d;
	    3*d  -4*d];
for i=1:N
  qam_sig(i,:)=mapping(dsource(i),:);
end
% received signal
for i=1:N
  n  = [sgma*randn sgma*randn];
  %[n(1) n(2)]=gngauss(sgma);
  r(i,:) = qam_sig(i,:) + n;
end
% detection and error probability calculation
numoferr=0;
for i=1:N
    
    % write your decoder here
    % The decision is named decis
    % Hint.
    % Given received signal vector r(i,:),
    % find the minimum distance among r(i,:) and the 16 signal points in
    % mapping(j,:), j=1,2,... 16. If the distance to mapping(k,:) is
    % minimum, decis=k
    min_distance = inf; % Initialize minimum distance to a large value
    for k = 1:16
        % Calculate Euclidean distance between received signal and each point in mapping
        distance = norm(r(i,:) - mapping(k,:));
        if distance < min_distance
            decis = k; % Update decision if current distance is smaller than previous minimum
            min_distance = distance;
        end
    end
  if (decis~=dsource(i))
    numoferr=numoferr+1;
  end
end
p=numoferr/(N);		
end