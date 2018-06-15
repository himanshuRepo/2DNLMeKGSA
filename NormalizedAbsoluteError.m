%===================================================================================
% MATLAB code for multi-level image thresholding segmentation using 2DNLMeKGSA.
% Author: Himanshu Mittal (himanshu.mittal224@gmail.com), 
%           Mukesh Saraswat (saraswatmukesh@gmail.com)
% Modified the code by:
%	Author : Athi Narayanan S
%	M.E, Embedded Systems,
%	K.S.R College of Engineering
%	Erode, Tamil Nadu, India.
%	http://sites.google.com/site/athisnarayanan/
%	s_athi1983@yahoo.co.in
%
% Developed in MATLAB R2015a
%
% Reference: "An optimum multi-level image thresholding segmentation using
%            non-local means 2D histogram and exponential Kbest gravitational 
%            search algorithm." Engineering Applications of Artificial 
%            Intelligence, Volume 71, Pages 226-235, Elsevier, 2018. 
%            https://doi.org/10.1016/j.engappai.2018.03.001
%
% File purpose: Measuring the Normalized Absolute Error Calculation.
%===================================================================================

function NAE = NormalizedAbsoluteError(origImg, distImg)

origImg = double(origImg);
distImg = double(distImg);

error = origImg - distImg;

NAE = sum(sum(abs(error))) / sum(sum(origImg));