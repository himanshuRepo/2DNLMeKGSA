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
% File purpose: "Main" file of the segmentation code.
%===================================================================================

clear;
close all;
clc;

runs=1;

method=[
'eKGSA';
];

for i=1:size(method,1)
    for j=1:runs
        methodname=strtok(method(i,:));
        filename=strcat(strtok(method(i,:)),'_',num2str(j),'run');
        
        [methodata]=callmain(methodname);
        
        
        xlswrite(filename,methodata);
        
        % For measuring the BDE, VoI, PRE, GCE
        seg_measure1(methodname,filename);
        % For measuring the other parameters, e.g. SSIM, RMSE, etc.
        seg_measure2(methodname,filename);
    end
end
