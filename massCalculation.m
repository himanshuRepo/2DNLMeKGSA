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
% File purpose: Computing the updated mass.
%===================================================================================

function [M]=massCalculation(fit,min_flag);

Fmax=max(fit); Fmin=min(fit); Fmean=mean(fit); 
[i N]=size(fit);

if Fmax==Fmin
   M=ones(N,1);
else
    
   if min_flag==1 %for minimization
      best=Fmin;worst=Fmax; %eq.17-18.
   else %for maximization
      best=Fmax;worst=Fmin; %eq.19-20.
   end
  
   M=(fit-worst)./(best-worst); %eq.15,

end

M=M./sum(M); %eq. 16.