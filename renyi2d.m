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
% File purpose: Objective function.
%===================================================================================

function [h]= renyi2d(xI,Pxy)
[M N]=size(xI);
alpha=0.5;

ind = Pxy == 0;
        ind = ind .* eps;
        Pxy = Pxy + ind;
        clear ind

    s1=xI(M,1);s2=xI(M,2);s3=xI(M,3);s4=xI(M,4);
    t1=xI(M,5);t2=xI(M,6);t3=xI(M,7);t4=xI(M,8);
    s1=round(s1);s2=round(s2);s3=round(s3);s4=round(s4);
    t1=round(t1);t2=round(t2);t3=round(t3);t4=round(t4);
    P0=0;
    for i=1:s1
        for j=1:t1
            P0=P0+Pxy(i,j);
        end
    end
    
    
    P1=0;
    for i=(s1+1):s2
        for j=(t1+1):t2
            P1=P1+Pxy(i,j);
        end
    end
    
    P2=0;
    for i=(s2+1):s3
        for j=(t2+1):t3
            P2=P2+Pxy(i,j);
        end
    end
    
        P3=0;
    for i=(s3+1):s4
        for j=(t3+1):t4
            P3=P3+Pxy(i,j);
        end
    end
    
        P4=0;
    for i=(s4+1):256
        for j=(t4+1):256
            P4=P4+Pxy(i,j);
        end
    end
        
    H0=0;
    for i=1:s1
        for j=1:t1
            sd=(Pxy(i,j)/P0)^alpha;
            H0 = H0+sd;                      
        end
    end
       
    H0=(1/(1-alpha))*log(H0+eps);
        
    H1=0;
    for i=(s1+1):s2
        for j=(t1+1):t2
            sd=(Pxy(i,j)/P1)^alpha;
            H1 = H1+sd;                      
        end
    end
        H1=(1/(1-alpha))*log(H1+eps);
    
    H2=0;
    for i=(s2+1):s3
        for j=(t2+1):t3
            sd=(Pxy(i,j)/P2)^alpha;
            H2 = H2+sd;                      
        end
    end
        H2=(1/(1-alpha))*log(H2+eps);
        
        
    H3=0;
    for i=(s3+1):s4
        for j=(t3+1):t4
            sd=(Pxy(i,j)/P3)^alpha;
            H3 = H3+sd;                      
        end
    end
        H3=(1/(1-alpha))*log(H3+eps);
    
    H4=0;
    for i=(s4+1):256
        for j=(t4+1):256
sd=(Pxy(i,j)/P4)^alpha;
            H4 = H4+sd;                      
        end
    end
    H4=(1/(1-alpha))*log(H4+eps);
    
h=H0+H1+H2+H3+H4;

   if ((isinf(h))|| (h==0))
        h=0;
    end
end