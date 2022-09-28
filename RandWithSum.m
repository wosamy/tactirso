function [ numbers ] = RandWithSum( targetSum, minVal, maxVal,n )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 
r = (maxVal-minVal).*rand(n,1) + minVal;
 
p= targetSum/sum(r);

% check=sum(r*p);
numbers=r*p;
end

