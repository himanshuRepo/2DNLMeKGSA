%===================================================================================
% MATLAB code for multi-level image thresholding segmentation using 2DNLMeKGSA.
% Author: Himanshu Mittal (himanshu.mittal224@gmail.com), 
%           Mukesh Saraswat (saraswatmukesh@gmail.com)
%
% Developed in MATLAB R2015a
%
% Reference: "An optimum multi-level image thresholding segmentation using
%            non-local means 2D histogram and exponential Kbest gravitational 
%            search algorithm." Engineering Applications of Artificial 
%            Intelligence, Volume 71, Pages 226-235, Elsevier, 2018. 
%            https://doi.org/10.1016/j.engappai.2018.03.001
%
% File purpose: Initializing the parameters of eKGSA algorithm.
%===================================================================================

function [gBest,gbestvalue,FEcount,etime,iteration] = maineKGSA(Lmax1,level,Pxy)


ElitistCheck=1; Rpower=1;
flag=0; % 1: minimization, 0: maximization

d = level-1;
NumberOfThresholdValues=2*d;
thdim=NumberOfThresholdValues;
Lmax=Lmax1;
XVmax = Lmax*ones(1,NumberOfThresholdValues);
XVmin = zeros(1,NumberOfThresholdValues);

Population_size=(thdim/2)*10;
max_it=1000;

%h      Defines the variant to be used
h=8;

objectiveFunction='renyi2d';
FE = 0;
tic

[Fbest,Lbest,FE,iteration]=eKGSA(objectiveFunction,Population_size,thdim,XVmin,XVmax,max_it,ElitistCheck,flag,Rpower,h,Pxy);
toc
gBest=Lbest;
gbestvalue=Fbest;
FEcount=FE;
etime = toc;
end