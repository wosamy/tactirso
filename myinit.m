%_________________________________________________________________________%
%  Rat Swarm Optimizer (RSO)                                              %
%                                                                         %
%  Developed in MATLAB R2019b                                             %
%                                                                         %

%_________________________________________________________________________%

function Pos=myinit(S,k,SearchAgents,dimension,upperbound,lowerbound)
Boundary= size(upperbound,2); 
if Boundary==1
    Pos=rand(SearchAgents,dimension).*(upperbound-lowerbound)+lowerbound;
end

if Boundary>1
    for i=1:dimension
        ub_i=upperbound(i);
        lb_i=lowerbound(i);
        Pos(:,i)=rand(SearchAgents,1).*(ub_i-lb_i)+lb_i;
    end
end



l=0;
for i=1:length(S) 
    if S(i).E>0
        l=l+1;
    end
end

p=k/l;
for r=1:ceil(.5*SearchAgents) % part  of  populTION say .5
    Pos(r,:)=zeros(1,dimension);
    for j=1:length(S)
        if (S(j).E>0)
            temp_rand=rand;
              Efactor=S(j).trust*(S(j).E/S(j).Eo)^2;
            v=(Efactor*(p/(1-p*mod(r,round(1/p)))));
            %Election of Cluster Heads
           if( temp_rand<= v)
                Pos(r,j)= (0.5).*rand(1,1)+0.5;%likely be ch
           else
               Pos(r,j)=(0.5).*rand(1,1);
            end
        end
    end
end
 


