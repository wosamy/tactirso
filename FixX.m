function [ newX ] = FixX( X,k )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
binX=GetBinForm(X);
chsIx=find(binX==1); % chs  has  1  value  and  0 for  cms
cmsIx=find(binX==0);
m=numel(chsIx);
if(m>k)
temp=randsample(chsIx,m-k);
binX(temp)=0;
end
if(m<k)
temp=randsample(cmsIx,k-m);
binX(temp)=1;
end
newX=binX;
end

