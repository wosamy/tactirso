%_________________________________________________________________________%
%  Rat Swarm Optimizer (RSO)                                              %
%                                                                         %
%  Developed in MATLAB R2019b                                             %
%                                                                         %
%  Designed and Developed: Dr. Gaurav Dhiman                              %
%                                                                         %
%         E-Mail: gdhiman0001@gmail.com                                   %
%                 gaurav.dhiman@ieee.org                                  %
%                                                                         %
%       Homepage: http://www.dhimangaurav.com                             %
%                                                                         %
%  Published paper: G. Dhiman et al.                                      %
%          A novel algorithm for global optimization: Rat Swarm Optimizer %
%               Jounral of Ambient Intelligence and Humanized Computing   %
%               DOI: https://doi.org/10.1007/s12652-020-02580-0           %
%                                                                         %
%_________________________________________________________________________%

function[Score,Position,Convergence]=rso2(S,k,Tthreshold,Resthreshold,Search_Agents,Max_iterations,Lower_bound,Upper_bound,dimension,objective)
Position=zeros(1,dimension);
Score=inf; 
Positions=myinit(S,k,Search_Agents,dimension,Upper_bound,Lower_bound);

Convergence=zeros(1,Max_iterations);
n=dimension;  
l=0;
x = 1;
y = 5;
R = floor((y-x).*rand(1,1) + x);
while l<Max_iterations
    for i=1:size(Positions,1)    
        Flag4Upper_bound=Positions(i,:)>Upper_bound;
        Flag4Lower_bound=Positions(i,:)<Lower_bound;
        Positions(i,:)=(Positions(i,:).*(~(Flag4Upper_bound+Flag4Lower_bound)))+Upper_bound.*Flag4Upper_bound+Lower_bound.*Flag4Lower_bound;
        newX=FixX(Positions(i,:),k);%randsample(k,1)
        fitness=objective(newX,S,Tthreshold,Resthreshold);
        
        %%%%NEW CODE FOR LOCAL SEARCH METHOD%%%%%
        rout_new=newX;
        randomnum=rand(1);
        if randomnum>0.5
            %%%%inversion algorithm
            ii1=randi(n);
            ii2=randi(n);
            rout_new(ii1:ii2)=rout_new(ii2:-1:ii1); 
            fitnessnew=objective(rout_new,S,Tthreshold,Resthreshold);
            %%%%end 
            if fitness>fitnessnew
                fitness=fitnessnew;
                newX=rout_new;
            end
            %%%%Single point exchange
            result1=randi(n);
            result2=randi(n);
            rout_new([result1 result2]) = rout_new([result2 result1]);
            fitnessnew=objective(rout_new,S,Tthreshold,Resthreshold);
            if fitness>fitnessnew
                fitness=fitnessnew;
                newX=rout_new;
            end
            %%%%End single point exchange
        end
        Positions(i,:)=newX;
        %%%%END NEW CODE LOCAL SEARCH METHOD%%%%%
        if fitness<Score 
            Score=fitness; 
            Position=Positions(i,:);
        end
        fitnessScores(i)=fitness;
    end
     
    A=R-l*((R)/Max_iterations); 
    
    for i=1:size(Positions,1)
        for j=1:size(Positions,2)     
            C=2*rand();          
            P_vec=A*Positions(i,j)+abs(C*((Position(j)-Positions(i,j))));                   
            P_final=Position(j)-P_vec;
            Positions(i,j)=P_final;
            
        end
    end
    l=l+1;    
    Convergence(l)=Score;
end
