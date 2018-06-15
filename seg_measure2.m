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
% File purpose: Measuring the SSIM, FSIM, MSE, RMSE, PSNR, NK, AD, SC, MD, 
%               and NAE values.
%===================================================================================

function seg_measure2(methodname,filename)

filePath = 'BW_dataset/';
directoryFiles = dir(filePath);
strc=strcat(methodname,'_result-set');
pathName = strc;
counter=1;
for fileIndex=1:length(directoryFiles)
    currentFilename = directoryFiles(fileIndex).name;
    disp(currentFilename);
    
    if (length(currentFilename)>4) && strcmp(currentFilename(end-3:end),'.mat')
        % It is an image file        
        imageName = currentFilename(1:end-4);
        dataFilename = [filePath imageName '.mat'];
        dataFilename1 = [pathName '/' imageName '.mat'];
        load(dataFilename);
        ref=refImage;
        load(dataFilename1);
        A=uint8(superpixelLabels);
        
        %% Calculate the SSIM
        [ssimval, ssimmap] = ssim(A,ref);
        SSIM(counter)=ssimval;

        %% Calculate the FSIM
        [fsim, fsimc] = FeatureSIM(ref, A);
        FSIM(counter)=fsim;

        %Mean Square Error 
MSE(counter) = MeanSquareError(ref, A);

        %Root Mean Square Error 
RMSE(counter) = sqrt(MeanSquareError(ref, A));

%Peak Signal to Noise Ratio 
PSNR(counter) = PeakSignaltoNoiseRatio(ref, A);

%Normalized Cross-Correlation 
NK(counter) = NormalizedCrossCorrelation(ref, A);

%Average Difference 
AD(counter) = AverageDifference(ref, A);

%Structural Content 
SC(counter) = StructuralContent(ref, A);

%Maximum Difference 
MD(counter) = MaximumDifference(ref, A);

%Normalized Absolute Error
NAE(counter) = NormalizedAbsoluteError(ref, A);

        counter=counter+1;
        
    end
end
A=[SSIM;FSIM;MSE;RMSE;PSNR;NK;AD;SC;MD;NAE];
st1=strcat(filename,'_II');
xlswrite(st1,A');
        %% Average
        AvgSSIM=mean(SSIM)        ;
        AvgFSIM=mean(FSIM);
        AvgMSE=mean(MSE);
        AvgRMSE=mean(RMSE);
        AvgPSNR=mean(PSNR);
        AvgNK=mean(NK);
        AvgAD=mean(AD);
        AvgSC=mean(SC);
        AvgMD=mean(MD);
        AvgNAE=mean(NAE);
B=[AvgSSIM;AvgFSIM;AvgMSE;AvgRMSE;AvgPSNR;AvgNK;AvgAD;AvgSC;AvgMD;AvgNAE];
st2=strcat(filename,'_avg_II');
xlswrite(st2,B);
end