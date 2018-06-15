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
% File purpose: Reading the image dataset, generating the 2D histogram and 
%               calling the eKGSA for optimal thresholds.
%===================================================================================

function [methodata]=callmain(methodname)

NORMALIZE_IMAGE = true;

method=methodname;
methodata=[];

filePath = 'Berkeley-dataset/';
directoryFiles = dir(filePath);
if strcmp(method,'eKGSA')
    pathName = 'eKGSA_result-set';
end
mkdir(pathName);



for fileIndex=1:length(directoryFiles)
    currentFilename = directoryFiles(fileIndex).name;
    disp(currentFilename);
    if (length(currentFilename)>4) && strcmp(currentFilename(end-3:end),'.jpg')
        % It is an image file
        imageName = currentFilename(1:end-4);
        dataFilename1 = [pathName '/' imageName '.mat'];
        
        imageFullPath = [filePath currentFilename];
        image = imread(imageFullPath);
        if ~isempty(whos('image'))
            [imageX,imageY,imageDepth] = size(image);
            if NORMALIZE_IMAGE
                % Resize the image to 320*240
                imageSize = 320;
                if imageX>=imageY
                    resizedImageY = round(imageY/imageX*imageSize);
                    image = imresize(image,[imageSize,resizedImageY]);
                else
                    resizedImageX = round(imageX/imageY*imageSize);
                    image = imresize(image,[resizedImageX,imageSize]);
                end
            end
            

            %% Convert to GreyScale
            image=rgb2gray(image);

            %% Defining the number of threshold levels
            partitionLevels=1; 
            if partitionLevels==1
                number_of_levels=5; 
            end
            
            %% Compute the non-local means (NLM) of an image
            I=image;
            level=number_of_levels;

            [m,n]=size(I);
            a=I;
            a0=im2double(a);
            t = 7;
            f = 2;

            a3=nlmeans(a0,t,f);

            for i=1:m
                for j=1:n
                    a4(i,j)=a3(i,j);
                end
            end

            a4 = (a4 - min(a4(:))) / (max(a4(:)) - min(a4(:)));
            a4 = im2uint8(a4);

            [m,n]=size(a);

            a0=(a);

            fxy=zeros(256,256);
            
            for i=1:m
                for j=1:n
                    c=a0(i,j);
                    d=(a4(i,j));
                    fxy(c+1,d+1)=fxy(c+1,d+1)+1;
                end
            end
            
            Pxy=fxy/m/n;
            %% Display the 2D-histogram     
            figure,
            mesh(Pxy);

            Lmax1=254;
            
            if size(I,3)==1 %grayscale image
                if strcmp(method,'eKGSA')
                    [gBest,gbestvalue,FEcount,etime,iteration]=maineKGSA(Lmax1,level,Pxy);
                end

                %% return optimal intensity
                intensity=round(gBest);     
                
                %% return fitness value
                fitness=gbestvalue;    
            end

            %% Extracting the threshold values
            Thresholds=intensity((number_of_levels):end);
            
            %% Sorting the thresholds
            srtThr=sort(Thresholds);

            %% Checking the case of similar threshold values
            v=find(diff(srtThr)==0);
            while (size(v,2) > 0)
                if size(I,3)==1 %grayscale image
                    if strcmp(method,'eKGSA')
                        [gBest,gbestvalue,FEcount,etime,iteration]=maineKGSA(Lmax1,level,Pxy);
                    end

                %% return optimal intensity
                intensity=round(gBest);     
                
                %% return fitness value
                fitness=gbestvalue;  

            end

            %% Extracting the threshold values
            Thresholds=intensity((number_of_levels):end);
            
            %% Sorting the thresholds
            srtThr=sort(Thresholds);

            %% Checking the case of similar threshold values
            v=find(diff(srtThr)==0);

            end
            
            X = imageGRAY(I,Thresholds);
            %% For displaying the segmented image
                         figure;
                         imshow(X);
            X1=X;
            X1=double(X1);
            
            
            %% for labelling
            X = imquantize(I,Thresholds);
            %% For displaying the segmented image
                         figure;
                         imshow(X);
            
            %% Post-Processing steps
            outdata={gBest(1),gBest(2),gBest(3),gBest(4),gBest(5),gBest(6),gBest(7),gBest(8),'' , gbestvalue, '' , FEcount, '' , etime, '' , iteration};
            methodata=[methodata; outdata];
            clearvars -except method X1 X J dataFilename1 filePath directoryFiles pathName fileIndex NORMALIZE_IMAGE methodata;
            sampleLabels=X;
            superpixelLabels=X1;
            save (dataFilename1, 'sampleLabels', 'superpixelLabels');
        end
    end
end

end
