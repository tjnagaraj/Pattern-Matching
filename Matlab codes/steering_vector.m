function [ a ] = steering_vector( pos )
%STEERING_VECTOR generates the steering vector for the given shifts and
%positon
%   @params   - "pos"   is the position
%             - "shifts" is the vector of shifts
%             - "a" is the steering vector output   D*1 vector

global shifts N;
a = exp((1*j*2*pi*pos*shifts)./N);
end

