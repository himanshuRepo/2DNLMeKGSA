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
% File purpose: Computing the updated velocity and positon.
%===================================================================================

function [X,V]=move(X,a,V)

%movement.
[N,dim]=size(X);
V=rand(N,dim).*V+a; %eq. 11.
V=fix(V);
X=X+V; %eq. 12.