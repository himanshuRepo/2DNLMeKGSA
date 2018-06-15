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
% File purpose: Computing the updated acceleration.
%===================================================================================

function a=eKGfield(M,X,G,Rnorm,Rpower,ElitistCheck,iteration,max_it,fit,min_flag,h);

[N,dim]=size(X);
final_per=2; %In the last iteration, only 2 percent of agents apply force to the others.


if ElitistCheck==1

    if h==8
           kbest=1.0055*((0.9944)^(iteration));
        kbest=round(N*kbest);
    
    end

else
    kbest=N; %eq.9.
end
[Ms ds]=sort(M,'descend');

for i=1:N
    E(i,:)=zeros(1,dim);
    for ii=1:kbest
        j=ds(ii);
        if j~=i
            R=norm(X(i,:)-X(j,:),Rnorm); %Euclidian distanse.
            for k=1:dim
                E(i,k)=E(i,k)+rand*(M(j))*((X(j,k)-X(i,k))/(R^Rpower+eps));
                %note that Mp(i)/Mi(i)=1
            end
        end
    end
end

%%acceleration
a=E.*G; %note that Mp(i)/Mi(i)=1


