%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%IRSO Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear;
clc;
flagfile=1;
efact=1;
%global Eelec;global Efs; global Emp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
seeds=[6060,2020,5050,4040,3030];
% seeds=(2020);
sink.xd=100;
sink.yd=100;
%Number of Nodes in the field
n=100;
x=0;
y=0;
%n=input('Enter the number of nodes in the space : ');
%Energy Model (all values in Joules)
%Initial Energy
Eo=0.5;
Eon=0.5;
Eoa=1;
MSize=4000;

%Eo=input('Enter the initial energy of sensor nJ : ');
%Field Dimensions - x and y maximum (in meters)
% xm=input('Enter x value for area plot : ');
% ym=input('Enter y value for area plot : ');
xm=100;
ym=100;

%x and y Coordinates of the Sink
% sink.xd=100;
% sink.yd=50;
Rc=40;
%Optimal Election Probability of a node
%to become cluster head
p=0.1;

%Eelec=Etx=Erx
Eelec=50*0.000000001;
%Transmit Amplifier types
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
%Data Aggregation Energy
EDA=5*0.000000001;
AR=0.50; % aggregation rate (AR)
%Computation of do
do=sqrt(Efs/Emp);

%Values for Hetereogeneity
%Percentage of nodes than are advanced or hetergenous
m=0.1;
%maximum number of rounds
%rmax=input('enter the number of iterations you want to run : ');
rmax=30000;


compromised=0.1;%%%%Percentage of compromised node
Tthreshold=0.5; %%%%Trust threshold
initialtotalenergy=(((n*m)*(Eoa)+(n-(n*m))*(Eon)));
% % Resthreshold=0.005; %%%REsidual energy threshold
% Resthreshold= (initialtotalenergy/initialtotalenergy)*(initialtotalenergy/n);
% %Resthreshold=(Eon*1)/100;


%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
ren=n*Eo-(m*n)*Eo;
if m==1
    ren=n*Eo;
end
randEn=RandWithSum(ren,0.2,0.8,m*n);

SHO_curveAll=[];

for run=1:1:5
    
    rng('default');
    rng(seeds(run)) %You can change the seed to have different positions for the node
    %     rng(20*run);
    countr=0;
    tottrustval=0;
    first_dead=[];
    DEAD=[];
    AvgEnergy=[];
    stdDev=[];
    half_dead=[];
    Fairness=[];
    PACKETS_TO_CH=[];
    PACKETS_TO_BS=[];
    STATISTICS=[];
    %Creation of the random Sensor Network
    ETotal=0;
    for i=1:1:n
        S(i).xd=rand(1,1)*xm;
        XR(i)=S(i).xd;
        S(i).yd=rand(1,1)*ym;
        YR(i)=S(i).yd;
        S(i).G=0;
        S(i).id=i;
        %initially there are no cluster heads only nodes 
        S(i).role='NA'; % node  role  CH or CM  intialy none NA 
        temp_rnd0=i;  
        %Random Election of Normal Nodes
        if (temp_rnd0>=m*n)
            S(i).E=Eo;
            S(i).Eo=Eo;
            S(i).MSize=MSize;
            S(i).type='N';% normal node
            
        end
        %Random Election of Advanced Nodes
        if (temp_rnd0<m*n+1)
            
            S(i).E=Eoa;%0.2 + (0.8-0.2).*rand(1,1); % [0.2,0.8]
            S(i).Eo= S(i).E;
            S(i).MSize=MSize;%randi([100,1900],1);%[100,1900]
            S(i).type='H';%heterogenous
            
            
        end
        S(i).di= sqrt( (S(i).xd-sink.xd )^2 + (S(i).yd-sink.yd )^2 );
        ETotal=ETotal+ S(i).Eo;
    end
    
    %%%RONDOMLY SET TRUST VALUES TO NODES%%%%%%%%%%
    %             compromised=0.1;%%%%Percentage of compromised node
    %             Tthreshold=0.5; %%%%Trust threshold
    trustset=randperm(n,(n*compromised));
    xmin=0.00001;
    xmax=Tthreshold;
    xfull=1;
    for th=trustset(1,:)
        S(th).trust=xmin+rand(1,1)*(xmax-xmin);
    end
    for i=1:n
        if S(i).trust>0
        else
            S(i).trust=xmax+rand(1,1)*(xfull-xmax);
        end
    end
    
    %%%%END ASSIGNING TRUST VALUE%%%%%%%%%%%%%%%%%%
    S(n+1).xd=sink.xd;
    S(n+1).yd=sink.yd;
    
    
    % flags  for  last  and  FND
    last_dead=0;
    flag_first_dead=0;
    flag_half_dead=0;
    TENPerRound=[];% total  engery  per round
    ema_en=[];
    ema_th=[];
    for r=1:1:rmax
        %counter for CHs
        countCHs=0;
        %counter for CHs per round
        rcountCHs=0;
        cluster=1;
        %Number of dead nodes
        dead=0;
        
        %counter for bit transmitted to Bases Station and to Cluster Heads
        packets_TO_BS=0;
        packets_TO_CH=0;
        %counter for bit transmitted to Bases Station and to Cluster Heads
        %per round
        PACKETS_TO_CH(r)=0;
        PACKETS_TO_BS(r)=0;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% get neighbor  and  distances
        for  i=1:1:n
            S(i).NB=[];
            S(i).disToNB=[];
            if S(i).E>0
                disToSink(i)=sqrt( (S(i).xd-(sink.xd) )^2 + (S(i).yd-(sink.yd) )^2 );
                S(i).disToSink=  disToSink(i);
                nj=1;
                neighbors=[];
                disToNeighbors=[];
                for j =1:1:n
                    if (i~=j && S(j).E>0)
                        dis=sqrt( (S(i).xd-(S(j).xd) )^2 + (S(i).yd-(S(j).yd) )^2 );
                        if(dis <Rc  )
                            neighbors(nj)=S(j).id;
                            disToNeighbors(nj)=dis;
                            nj=nj+1;
                        end
                    end
                end
                
                S(i).NB=neighbors;
                S(i).disToNB= disToNeighbors;
            else % dead  node
                S(i).E=0;
            end
        end
        %%%%%%%%%%%%%%%%%% sink neigbors
        sink.NB=[];
        for  i=1:1:n
            if(S(i).disToSink <Rc )
                sink.NB=union (sink.NB,S(i).id);
            end
        end
        for i=1:1:n
            S(i).no_of_CMs=0;
            %checking if there is a dead node
            if (S(i).E<=0)
                dead=dead+1;
                S(i).E=0;S(i).trust=0;
            end    
        end
        % when  all  node  dies
        if dead==n
            fprintf(' last  node  die at:%d \n',r);
            last_dead=r;
            break;
        end 
        STATISTICS(r).DEAD=dead;
        DEAD(r)=dead;
        %When the first node dies
        if (dead>=1)
            if(flag_first_dead==0)
                first_dead=r;
                flag_first_dead=1;
            end
        end
        %When the half node dies
        if (dead>=n/2)
            if(flag_half_dead==0)
                halfn_dead=r;
                flag_half_dead=1;
            end
        end
        if (dead>(n/2)-5 && dead<(n/2)+5)
            half_dead= union(half_dead, r);
        end 
        for i=1:1:n
            S(i).role='NA';
        end 
        % election  process   occurs  at  first round  or  when  total network energy below threshold value
        
        %%%%%%%%%%%%%%%%%%%%%%  Cluster_head_election%%%%%%%%%%%%%%%%%%%%%%%5
        liveNodes=0;totalLiveEnergy=0;
        for i=1:1:n
            if (S(i).E>0)
                liveNodes=liveNodes+1;
                totalLiveEnergy=totalLiveEnergy+S(i).E;
            end
            
        end
        TENPerRound(r)=totalLiveEnergy/(liveNodes);
        ema_en(r)=TENPerRound(r);
        ema_th(r)=(totalLiveEnergy/ETotal)*Tthreshold;%  change  rate  in egergy * trust threshold
     
        Resthreshold=ema_en(r);%TENPerRound(r);%totalLiveEnergy/(2*liveNodes);
        Tthreshold_withENrate=ema_th(r);

        if r>10
        %%%https://www.investopedia.com/terms/e/ema.asp 
        %for EMA EMA = Closing price x multiplier + EMA (previous day) x (1-multiplier)
        %  multiplier=[2 ÷ (number of observations + 1)]
        pastr=10;
        multiplier =2/(pastr+1);
        ema_en(r)=TENPerRound(r)*multiplier+ema_en(r-1)*(1-multiplier);
        ema_th(r)=Tthreshold_withENrate*multiplier+ema_th(r-1)*(1-multiplier);
        Resthreshold=ema_en(r);
         Tthreshold_withENrate=ema_th(r);
        end
        % Load details of the selected benchmark function
        k=ceil(p*liveNodes);%  number of reqired chs and  uodate  it  dynamically  according  to  number  of  live  nodes

        SearchAgents=30;
        Max_iterations=150;
        lowerbound=0;
        upperbound=1;
        dimension=n;
        fitness=@Fit2;
        [Best_score,Best_pos,SHO_curve]=rso2(S,k,Tthreshold_withENrate,Resthreshold,SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness);

        SHO_curveAll= [SHO_curveAll;SHO_curve];       
        C=[];
        tottrust=0;
        CHSARE=find(GetBinForm(Best_pos)==1);
        countCHs=numel(CHSARE);
        fprintf('round= %d - ch count=%d lives=%d\n',r,countCHs, liveNodes);
        if countCHs>0
            for io=CHSARE(1,:)
                S(io).role='CH';
                tottrust=S(io).trust+tottrust;
                S(io).G=round(1/p)-1;
                C(cluster).xd=S(io).xd;
                C(cluster).yd=S(io).yd;
                distance=sqrt( (S(io).xd-(S(n+1).xd) )^2 + (S(io).yd-(S(n+1).yd) )^2 );
                C(cluster).distance=distance;
                C(cluster).id=io;
                C(cluster).no_of_CMs=0;
                C(cluster). CMsBits= S(io).MSize;
                C(cluster).TotalCMsBits= S(io).MSize;
                cluster=cluster+1;
            end
        end
        
        %%%%RAT END CH SELECTION%%%%%%%%%%%%%%%%
        
        numberofcluster=cluster-1;
        if  numberofcluster>0
            averagetust=(tottrust/numberofcluster);
        else
            averagetust=0;
        end
        alivecount=0;
        for i=1:n
            if S(i).E>0
                alivecount=alivecount+1;
            end
        end  
        for i=1:1:n
            n1data=[];
            if ( strcmp(S(i).role,'NA') && S(i).E>0 )
                    if(cluster-1>=1)
                        min_dis=sqrt( (S(i).xd-(S((C(1).id)).xd) )^2 + (S(i).yd-(S((C(1).id)).yd) )^2 );
                        min_dis_cluster=1;
                        for c=1:1:cluster-1% select nearest  CH to the node
                            temp=sqrt( (S(i).xd-(S((C(c).id)).xd) )^2 + (S(i).yd-(S((C(c).id)).yd) )^2 );
                            if(temp< min_dis)
                                min_dis=temp;
                                min_dis_cluster=c;
                            end
                        end              
                        S(i).role='CM';
                        cm= S(C(min_dis_cluster).id).no_of_CMs;
                        S(C(min_dis_cluster).id).no_of_CMs=cm+1;
                        C(min_dis_cluster).no_of_CMs=C(min_dis_cluster).no_of_CMs+1;
                        C(min_dis_cluster). CMsBits=union(C(min_dis_cluster).CMsBits,S(i).MSize);
                        C(min_dis_cluster).TotalCMsBits=C(min_dis_cluster).TotalCMsBits+S(i).MSize;
                        %Energy dissipated by associated Cluster Head
                        min_dis=sqrt( (S(i).xd-(S((C(min_dis_cluster).id)).xd) )^2 + (S(i).yd-(S((C(min_dis_cluster).id)).yd) )^2 );
                    
                        MSize=S(i).MSize;
                        if (min_dis>do)
                            S(i).E=S(i).E- ( Eelec*(MSize) + Emp*MSize*( min_dis * min_dis * min_dis * min_dis));
                        end
                        if (min_dis<=do)
                            S(i).E=S(i).E- ( Eelec*(MSize) + Efs*MSize*( min_dis * min_dis));
                        end
                    
                        %Energy dissipated of CH for receiving
                        S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (Eelec + EDA)*MSize );
                        packets_TO_CH=packets_TO_CH+1;                 
                    
                    else % no  head  for  this  round send  direct  to  BS
                        MSize=S(i).MSize;
                    
                        min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
                        if (min_dis>do)
                            S(i).E=S(i).E- ( Eelec*(MSize) + Emp*MSize*( min_dis * min_dis * min_dis * min_dis));
                        end
                        if (min_dis<=do)
                            S(i).E=S(i).E- ( Eelec*(MSize) + Efs*MSize*( min_dis * min_dis));
                        end
                        packets_TO_BS=packets_TO_BS+1;
                    end
            end
        end       
        for c=1:1:cluster-1
            %Calculation of Energy dissipated for CHs tosend to BS          
            i=C(c).id;           
            MSize=S(i).MSize;           
            distance=C(c).distance; 
            if (distance>do)
                S(i).E=S(i).E- ( (Eelec)*(MSize) + Emp*MSize*( distance*distance*distance*distance ));
            end
            if (distance<=do)
                S(i).E=S(i).E- ( (Eelec)*(MSize)  + Efs*MSize*( distance * distance ));
            end
            packets_TO_BS=packets_TO_BS+1;  
        end
        %%% statistic        
        total=0.0;
        Energy=[];EnergyRatio=[];
        count=0;
        for i=1:1:n
            if(S(i).E>0)
                Energy(i)=S(i).E;EnergyRatio(i)=S(i).E/S(i).Eo;
                total=total + Energy(i);
                count=count+1;
            else
                Energy(i)=0;EnergyRatio(i)=0;
            end
        end
        average =total/n;
        averageremenergy=total/n;
        numberofalivenodes=count;        
        TEnergy(r)=total;% total remaining
        AveEnergy(r)=average;
        stdDev(r)=std( Energy);
        EnMap(r,:)= Energy';
        STATISTICS(r).RemainingEnergy=Energy;
        STATISTICS(r).EnergyRatio=EnergyRatio;
        STATISTICS(r).EDR=abs(max(EnergyRatio)-min(EnergyRatio(EnergyRatio>0))); 
        totalCMs=0;  
        PACKETS_TO_CH(r)=packets_TO_CH;
        PACKETS_TO_BS(r)=packets_TO_BS;
        STATISTICS(r).CLUSTERHEADS=cluster-1;
        STATISTICS(r).CMs=totalCMs;
         if flagfile==1
            filename=     strcat('Runs','_rat_',datestr(now,'mm-dd-yyyy HH-MM'),'.txt');
            fileID = fopen(filename, 'w');
                         fprintf(fileID,'compromised= %f ,Sink ( %d,%d), Eo=%f  , Eoa= %f ,Tthreshold=%f\n\n',compromised,sink.xd,sink.xd,Eo,Eoa,Tthreshold); 

            fprintf(fileID,'Run,Round,numberofcluster,averagetust, averageremenergy,numberofalivenodes,TotalEn\n\n'); 
            flagfile=0;
        end
        fileID = fopen(filename, 'a');
        fprintf(fileID,'%d, %d ,  %d , %6.4f, %12.8f , %d , %6.4f \n',run,r,numberofcluster,averagetust, averageremenergy,numberofalivenodes, total);
        fclose(fileID);
        if averagetust>=0
            countr=countr+1;
            tottrustval=tottrustval+averagetust;
        end     
    end
    trustfinal=tottrustval/countr;     
    fileIDSummary = fopen(strcat('Summary ','_rat_',datestr(now,'mm-dd-yyyy HH'),'.txt'), 'a');
    fprintf(fileIDSummary,'compromised= %f ,Sink ( %d,%d), Eo=%f  , Eoa= %f ,Tthreshold=%f\n\n',compromised,sink.xd,sink.xd,Eo,Eoa,Tthreshold); 
    fprintf(fileIDSummary,'Run=\t %d \t FND=\t %d \t HND= \t %d \t LND=\t  %d  \t,  Trust=\t %f \tfor  current run and new  run  started\n',run,first_dead,halfn_dead, last_dead, trustfinal);
    fclose(fileIDSummary);
    fprintf('First  node  die:\t %d \t HND= \t %d \t and last die\t  %d  \t,  Trust=\t %f \tfor  current run and new  run  started\n',first_dead,halfn_dead, last_dead, trustfinal);
    
end
