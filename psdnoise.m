%===================================================================================
% MATLAB code for multi-level image thresholding segmentation using 2DNLMeKGSA.
% Author: Himanshu Mittal (himanshu.mittal224@gmail.com), 
%           Mukesh Saraswat (saraswatmukesh@gmail.com)
% Modified this file for the non-commercial purpose only.
% (https://quasar.ugent.be/bgoossen/download_nlmeans/)
%
% Developed in MATLAB R2015a
%
% Reference: "An optimum multi-level image thresholding segmentation using
%            non-local means 2D histogram and exponential Kbest gravitational 
%            search algorithm." Engineering Applications of Artificial 
%            Intelligence, Volume 71, Pages 226-235, Elsevier, 2018. 
%            https://doi.org/10.1016/j.engappai.2018.03.001
%
% File purpose: Part of the file that computes the non-local means of an image.
%===================================================================================

function [Noise, PS] = psdnoise(M, N, type, s)
%PSDNOISE: Generates noise using power spectral density function
%   (C) 2006 Bart Goossens, Ghent University
%   Department of Telecommunications and Information Processing
%   UGent-TELIN-IPI-IBBT
%   bart.goossens (at) telin.ugent.be
% Arguments:
% - M:     height of the image
% - N:     width of the image
% - type:  noise preset type (each type represents a different power spectral
%          density - PSD)
% - s:     spatial scaling factor for the noise (optional, default: s=1)
% Outputs:
% - Noise: the generated noise
% - PS:    the power spectral density function of the noise
if nargin == 0
    standalone = 1;
    M = 256; N = 256; type=1;
else
    standalone = 0;
end;

if nargin < 4
   s = 1;
end;    

% actually, this can also be implemented using the meshgrid function
u = ((-M/2):(M/2-1))'*ones(1,N)/M/s;
v = ones(M,1)*((-N/2):(N/2-1))/N/s;

switch(type)
    case 0
        P = ones(M,N);
    case 1
        P = 0.05+exp(-30* ((u-0.02).^2)) .* exp(-3* ((v-0.02).^2));
    case 2
        P = 0.005 + 4*exp(-30* ((u-0.17).^2 + (v-0.10).^2));
    case 3
        % circular symmetrical psd; found on digital pictures
        P = 0.005+exp(-10* (u.^2 + v.^2));
    case 4
        P = 0.05+exp(-300* ((u-0.02).^2)) .* exp(-30* ((v-0.02).^2));
    case 5
        P = 0.005 + exp(-60* ((u-0.17).^2 + (v-0.10).^2));
    case 6
        P = 0.05+ (exp(-20 *((u - v).^2 )));
    case 7
        P = 0.05 + (sqrt(u.^2 + v.^2) < 0.1);
    case 8
        P = 0.05 + (((v).^2 + (0.1*u+7*v).^2) < 0.1);
    case 9
        % ideal bandpass noise
        theta = -3*pi/8; % rotation angle
        u2 = cos(theta) * u - sin(theta) * v;
        P = 0.05 + (abs(u2) > 0.005) .* (abs(u2) < 0.06);
    case 10
        % bandpass noise
        theta = 3*pi/8; % rotation angle
        u2 = cos(theta) * u - sin(theta) * v;
        v2 = sin(theta) * u + cos(theta) * v;
        P = 0.05 + exp(-abs(v2.^2+u2.^2-0.05)*60);    
    case 11
        P = 0.05 + exp(-(sqrt(u.^2+v.^2) - 0.05).^2*8000);
    case 12
        P = 0.05 + abs(real(exp(-3000*((u-0.17).^2+(v-0.1).^2))));
    case 13
        % line pattern noise; found on analogue video
        P = exp(-4000 * ((u-0.1).^2 + (v-0.12).^2));
    case 13.1
        % line pattern noise (modified); found on analogue video
        P = 1e-3+exp(-4000 * ((u-0.1).^2 + (v-0.12).^2));
    case 13.2
        % line pattern noise (modified); found on analogue video
        P = 1e-3+exp(-4000 * ((u+0.1).^2 + (v-0.12).^2));
        case 14
        % line pattern noise; found on analogue video
        P = 1-exp(-0.1* (u.^2 + v.^2))  +1e3* exp(-4000 * ((u-0.1).^2 + (v-0.12).^2));               
    case 15
        % double line pattern noise; found on analogue video
        P = exp(-2000 * ((u-0.1).^2 + (v-0.12).^2)) + ...
            exp(-3000 * ((u+0.15).^2 + (v-0.22).^2)) + 1e-3;   
 end;        

% symmetrical extension (to make the filter 'real'-valued)
P = fftshift(sqrt(P));
P = 0.5 * (P + circshift(P(end:-1:1,end:-1:1),[1 1])); % trick :)

if standalone
    figure;
    imagesc(abs(fftshift(P)));colormap(gray);
    title('Power spectral density (PSD)');
end;    

x = fft2(randn(M, N));
x = x ./ abs(x);
x = x .* (P);

% normalization
P = P / sqrt(mean2(P.^2));

Noise = ifft2(x);
Noise = Noise / sqrt(mean2(Noise.^2));


if nargout>1    
    PS = P.^2;
end;    
