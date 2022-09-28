function [ binX ] = GetBinForm( X )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
  
            
                for j=1:numel(X)
%                     sig=( exp( (X(j)))-1)/( exp( (X(j)))+1)
                    if X(j)>0.5
                        binX(j)=1;
                    else
                        binX(j)=0;
                    end
                end


end

