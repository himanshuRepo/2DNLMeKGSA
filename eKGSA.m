%===================================================================================
% MATLAB code for multi-level image thresholding segmentation using 2DNLMeKGSA.
% Author: Himanshu Mittal (himanshu.mittal224@gmail.com), 
%           Mukesh Saraswat (saraswatmukesh@gmail.com)
% Modified the gravitational search algorithm code by Esmat Rashedi, 2010.
%
% Developed in MATLAB R2015a
%
% Reference: "An optimum multi-level image thresholding segmentation using
%            non-local means 2D histogram and exponential Kbest gravitational 
%            search algorithm." Engineering Applications of Artificial 
%            Intelligence, Volume 71, Pages 226-235, Elsevier, 2018. 
%            https://doi.org/10.1016/j.engappai.2018.03.001
%
% File purpose: Calling the eKGSA algorithm for optimal thresholds.
%===================================================================================

function [Fbest,Lbest,FE,iteration]=eKGSA(fname,N,thdim,low,up,max_it,ElitistCheck,min_flag,Rpower,h,Pxy)

%V:   Velocity.
%a:   Acceleration.
%M:   Mass.  Ma=Mp=Mi=M;
%dim: Dimension of the test function.
%N:   Number of agents.
%X:   Position of agents. dim-by-N matrix.
%R:   Distance between agents in search space.
%[low-up]: Allowable range for search space.
%Rnorm:  Norm in eq.8.
%Rpower: Power of R in eq.7.


 Rnorm=2; 
 dim=thdim/2;
 FE=0;
 FEmax=dim*1000;


%% Random initialization for agents.

dim1=dim;
X1=initialization(dim1,N,up,low); 
for si=1:length(X1)
           X1(si,:)=sort(X1(si,:)); 
end
dim2=dim;
X2=initialization(dim2,N,up,low); 
for si=1:length(X2)
           X2(si,:)=sort(X2(si,:)); 
end


X=[X1 X2];


%create the best so far chart and average fitnesses chart.
BestChart=[];MeanChart=[];

V=zeros(N,thdim);

for iteration=1:max_it
%     iteration
    
    %Checking allowable range. 
    
    X=space_bound(X,up,low);
    
    Xtemp1=X(:,1:dim);
for si=1:length(Xtemp1)
           Xtemp1(si,:)=sort(Xtemp1(si,:)); 
end

Xtemp2=X(:,(dim+1):thdim);
for si=1:length(Xtemp2)
           Xtemp2(si,:)=sort(Xtemp2(si,:)); 
end
    
   X=[Xtemp1 Xtemp2] ;    
    

    for i=1:N
    [fitness(i)]=feval(fname,X(i,:),Pxy);
    FE=FE+1;
    end
    
    if min_flag==1
    [best best_X]=min(fitness); %minimization.
    else
    [best best_X]=max(fitness); %maximization.
    end        
    
    if iteration==1
       Fbest=best;Lbest=X(best_X,:);
    end
    if min_flag==1
      if best<Fbest  %minimization.
       Fbest=best;Lbest=X(best_X,:);
      end
    else 
      if best>Fbest  %maximization
       Fbest=best;Lbest=X(best_X,:);
      end
    end
      

%Calculation of M. eq.14-20
[M]=massCalculation(fitness,min_flag); 

%Calculation of Gravitational constant. eq.13.
G=Gconstant(iteration,max_it); 

%Calculation of accelaration in gravitational field. eq.7-10,21.
a=eKGfield(M,X,G,Rnorm,Rpower,ElitistCheck,iteration,max_it,fitness,min_flag,h);

%Agent movement. eq.11-12
[X,V]=move(X,a,V);

end %iteration

end

