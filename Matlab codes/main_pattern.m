%% Program to evaluate the performance of Sparse Graph Code based search algorithm

clc
close all
clear all

%% Parameters

global N_db N eta M N1 N2 N3 N4 N5 N6 N7 N8 P num_stages num_shifts
k=1;
eta = 0;

[signal, Fs1] = audioread('speech.wav');

N = 4.8*10^6

for l = 1
    l
start =tic;


% N_tilde = 10^8;
num_blocks = 1;
N_tilde = N;

%M = 10^5;
num_locs = 1;
%M = length(y);

%K = floor(eta*M)

% RSIDFT parameters 
num_stages =2;


if(l ==1)
N1=2150;
N2=2157;
N3=1;
N4=1;
N5=1;
N6=1;
N7=1;
N8=1;
end
P = floor((N_tilde)/(N1*N2*N3*N4*N5*N6*N7*N8));

N= P*N1*N2*N3*N4*N5*N6*N7*N8;
M = 5*Fs1;
x = signal(1:N,1);       % 100s clip
y = x(1*Fs1:6*Fs1-1,1);    % music clip
%y_edit = [y(1:floor(length(y)/2),1) y(floor(length(y)/2)+1:length(y),1)];

len_y = length(y)
%len_y_edit = length(y_edit)

N_db = N * num_blocks;


alpha = log10((N_db)/(N1))/log10(N_db)
alpha_block = log10((N)/(N1))/log10(N);
mu = log10(M)/log10(N_db)

num_shifts =10;
% num_shifts = ceil(log(5*N))

scheme=0;               % scheme : '0' -  less-sparse (alpha = 1/3)     '1' -   very-sparse (alpha=2/3)

%% Test Program

%num_shifts = 10

samples_gain = N_db/((num_blocks * num_stages * num_shifts*((N)^alpha_block)))
% computational_gain = (N_db * log(N_db))/  


num_trials=1;
Xs_match=[];
Xs_mismatch=[];
num_checks=0;

prob_failure=0;

locs = loc_generation(N_db,M,num_locs);
%query = sign(randi([0,1],1,M)-0.5);
query = y;

error.locs{l,1}= alpha;
error.locs{l,2} =locs;
avg_checks_vis = 0;

for i =1:num_trials
    if(mod(i,10)==0)
 fprintf('Block number %d, time elapsed %d sec, perror %d \n',i, toc(start), prob_failure/(i-1));
    end
  loc_i = locs(find(((i-1)*N < locs) &(locs < (i*N - M))))-((i-1)*N);

start1 = tic;
start2 = tic;

n_locs_i = numel(loc_i);

input = fft(x,N) .* conj(fft(query,N));

figure;

plot(ifft(input))

[X_est, Xs, num_vis_checks, num_checks]= FFAST_robust_3stages( scheme, input);

 
loc_est = find(X_est~=0)

if(n_locs_i ~=0)
num_missed =n_locs_i - (numel(intersect(loc_est, loc_i )));
if length(loc_est) > (n_locs_i - num_missed) 
    error.events{k,1} = alpha;
    error.events{k,2} = loc_i;
    error.events{k,3} = loc_est;
    error.events{k,4} = 'wrong';
    k=k+1;
elseif(num_missed >0)
    error.events{k,1} = alpha;
    error.events{k,2} = loc_i;
    error.events{k,3} = loc_est;
    error.events{k,4} = 'missed';
    k=k+1;
end
prob_failure=(num_missed/n_locs_i)+ prob_failure;

end
avg_checks_vis = avg_checks_vis + num_vis_checks;
end
perror(l) = prob_failure/num_trials;
avg_checks_vis(l) = avg_checks_vis/num_trials;

error.locs{l,4} = perror(l); 
error.locs{l,3} = [avg_checks_vis(l) num_checks]
end
       
        