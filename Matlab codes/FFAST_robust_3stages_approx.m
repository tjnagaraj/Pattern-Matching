%FFAST Peeling decoder with 3-stages
function [ X_est, Xs_bin, num_vis_checks, num_checks ] = FFAST_robust_3stages(scheme, x )

%% Parameters

% N = N1* N2 * N3,  where N1, N2, N3 are co-prime
% scheme : '0' means less-sparse and '1' means very-sparse
% gamma - parameter choosen to adapt threshold value to noise
% x - input signal


%%  Params for the simulation
global P N N1 N2 N3 N4 N5 N6 N7 N8 eta Cconn shifts num_shifts num_stages

N_i = [N1 N2 N3 N4 N5 N6 N7 N8];


%% --- Robust FFAST implementation
xs=[];
Xs=[];


% for i=1:length(shifts)
%     xs(i,:) = circshift(x,[0 -1*shifts(1,i)]);       % shifted versions of x
% end

if(scheme==0)  % less-sparse regime
    f=N./N_i;    % # of samples in the each stage per branch
else
    f=P*N_i;     % # of samples in the each stage per branch
end

d= N./f;       %down-sampling factors

shifts = [0 randperm(N,num_shifts-1)];
xs_down= {};

num_checks = sum(f(1:num_stages));

for i=1:num_stages
    shift_mat = [];
    sampling_mat=[];
    shift_mat = repmat(transpose(-1* shifts),[1,f(i)]);
    sampling_mat = mod(shift_mat +  repmat(0:d(i):(f(i)-1)*d(i),[num_shifts,1]),N);
    xs_down{i,1} = x(sampling_mat+1);
end

Xs ={};
Xs_bin={};


for i=1:num_stages
    Xs{i,1} = ifft(xs_down{i,1},f(i),2);
    Xs_bin{i,1} = transpose(Xs{i,1});
end

% Graph Connections
Cconn={};
for stag_num = 1:num_stages
    for bin_num = 1:f(stag_num)
        Cconn{stag_num,bin_num}= bin_num-1:f(stag_num):bin_num-1+(d(stag_num)-1)*f(stag_num);
    end
end

%% Peeling process

X_est =zeros(1,N) ;     % Estimated signal in Fourier domain
l=50;                     % Number of iterations required to recover the signal

op= -100*ones(num_stages,max(f));
value= -100*ones(num_stages,max(f));
pos= -100*ones(num_stages,max(f));

num_vis_checks=0;


for stag_num = 0:num_stages-1
    for bin_num = 0:f(stag_num+1)-1
        y = transpose(Xs_bin{stag_num+1,1}(bin_num+1,:));
        if i==1
            [op1, value1 , pos1, pos_flag] = bin_detector(y,stag_num,bin_num);
            op(stag_num+1,bin_num+1)=op1;
            value(stag_num+1,bin_num+1) = value1;
            pos(stag_num+1,bin_num+1)=pos1;
            if(pos_flag==1)
                num_vis_checks = num_vis_checks+1;
            end
        end
    end
end



for i=2:l
    %  iter=i
    sing_flag =0;
    for stag_num = 0:num_stages-1
        for bin_num = 0:f(stag_num+1)-1
            
            if(op(stag_num+1,bin_num+1)==0)
                %                 fprintf('%s,%d,%s,%d,%s\n','The checknode(',stag_num,',',bin_num, ') is a Zero-ton' )
            elseif(op(stag_num+1,bin_num+1)==1)
                sing_flag=1;
                %      fprintf('%s %d %s %d %s\n','The checknode(',stag_num,',',bin_num, ') is a Single-ton' )
                X_est(pos(stag_num+1,bin_num+1)+1)=value(stag_num+1,bin_num+1);
                for snum = 1: num_stages
                    Cind_init = mod(pos(stag_num+1,bin_num+1),f(snum));
                    
                    if(eta~=0)
                        Cind=[];
                        for cind = Cind_init
                            if(op(snum,cind+1)~=3)
                                Cind = [Cind cind];
                            end
                        end
                    else
                        Cind = Cind_int;
                    end
                    
                    Xs_bin{snum,1}(Cind+1,:)= Xs_bin{snum,1}(Cind+1,:)-value(stag_num+1,bin_num+1)*steering_vector(pos(stag_num+1,bin_num+1));
                    y = transpose(Xs_bin{snum,1}(Cind+1,:));
                    [op1, value1 , pos1, pos_flag] =bin_detector(y,snum-1,Cind);
                    op(snum,Cind+1) = op1;
                    value(snum,Cind+1) = value1;
                    pos(snum,Cind+1) = pos1;
                    if pos_flag==1
                        num_vis_checks = num_vis_checks+1;
                    end
                    
                end
            else
                %                 fprintf('%s,%d,%s,%d,%s\n','The checknode(',stag_num,',',bin_num, ') is a Multi-ton')
            end
            
        end
    end
    if(sing_flag ==0)
        break;
    end
end

end

