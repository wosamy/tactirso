% F20

function v = Fit2(X,S,Tthreshold,Resthreshold)

binX=GetBinForm(X);
heads=find(binX==1);

if isempty(heads)
    v=inf;
else
    v=0;vt=0;ve=0;vn=0;
    
    %          fln=0;
    flrt=0; % for energy
    fltt=0;% for trust
    fdi=0;
    % chech if there is chs separation
%     for j=1:length(heads)
%         for k=1:length(heads)
%             if j~=k
%                 di= sqrt( (S(j).xd-S(k).xd )^2 + (S(j).yd-S(k).yd )^2 );
%                 if(di>40)
%                     fdi=fdi+1;% how many chs sperated
%                 end
%             end
%         end
%     end
          
    
    for j=heads(1,:)
         ve=ve+ (1/(1+S(j).E));vt=vt+(1/(1+S(j).trust));
        if S(j).E>=Resthreshold
            flrt=flrt+1;
            
        end
        
        if Tthreshold<=S(j).trust
            fltt=fltt+1;
        end
    end
    
    v=(ve*(1/(flrt+1))+(1/(1+fltt))*vt);%+(1/(1+fdi));
      %v=((1/(flrt+1))+(1/(1+fltt)));%+(1/(1+fdi));
    
    
    
    
    
    
end


end

