%===================================================================================
% MATLAB code for multi-level image thresholding segmentation using 2DNLMeKGSA.
% Author: Himanshu Mittal (himanshu.mittal224@gmail.com), 
%           Mukesh Saraswat (saraswatmukesh@gmail.com)
% Modified this file for the non-commercial purpose only.
% (https://people.eecs.berkeley.edu/~yang/software/lossy_segmentation/)
%
% Developed in MATLAB R2015a
%
% Reference: "An optimum multi-level image thresholding segmentation using
%            non-local means 2D histogram and exponential Kbest gravitational 
%            search algorithm." Engineering Applications of Artificial 
%            Intelligence, Volume 71, Pages 226-235, Elsevier, 2018. 
%            https://doi.org/10.1016/j.engappai.2018.03.001
%
% File purpose: Measuring the BDE, VoI, PRE, and GCE values.
%===================================================================================
 
function seg_measure1(methodname,filename)
benchPath = 'benchset/';
strc=strcat(methodname,'_result-set/')
testPath = strc;
cropBenchImage = true;

testImageList = dir(testPath);
if isempty(testImageList)
    error('Cannot find the image directories.');
end

resizeScale = 320;
averageBoundaryError = 0;
averageRI = 0;
averageVOI = 0;
averageGCE = 0;

imageCount = 0;
for imageIndex=1:length(testImageList)
    if testImageList(imageIndex).name(1)=='.'
        continue;
    end
    
    testFilename = [testPath testImageList(imageIndex).name];
    if ~strcmp(lower(testFilename(end-3:end)),'.mat')
        % Not a valid data file
        continue;
    end
    
    benchFilename = [benchPath testImageList(imageIndex).name];
    if isempty(dir(benchFilename))
        % The bench data for this image is not available
        warning(['Cannot find the bench file for ' testImageList(imageIndex).name]);
    end
    
    disp(['Processing ' testImageList(imageIndex).name]);
    superpixelLabels = [];
    load(benchFilename);
    load(testFilename);
    imageCount = imageCount + 1;
    
    % Comparison script
    totalBoundaryError = 0;
    sumRI = 0;
    sumVOI = 0;
    sumGCE = 0;

    [imageX, imageY] = size(sampleLabels);
    [benchX, benchY] = size(imageLabelCell{1});
    
    for benchIndex=1:length(imageLabelCell)
        benchLabels = imageLabelCell{benchIndex};
        
        % Rescale benchimage
        if benchX>benchY
            benchY=benchY*resizeScale/benchX;
            benchX=resizeScale;
        else
            benchX=benchX*resizeScale/benchY;
            benchY=resizeScale;
        end
        benchLabels=imresize(benchLabels, [benchX, benchY],'nearest');
        
        cropBoundarySize = (size(benchLabels,1)-size(sampleLabels,1))/2;
            
        % update the four error measures:        
        totalBoundaryError = totalBoundaryError + compare_image_boundary_error(benchLabels, sampleLabels);        
        
        [curRI,curGCE,curVOI] = compare_segmentations(sampleLabels,benchLabels);       
        sumRI = sumRI + curRI;
        sumVOI = sumVOI + curVOI;
        sumGCE = sumGCE + curGCE;        
    end
    BE(imageCount)=totalBoundaryError;
    RI(imageCount)=sumRI;
    VOI(imageCount)=sumVOI;
    GCE(imageCount)=sumGCE;
    % update the averages... note that sumRI / length(imageLabelCell) is
    % equivalent to the PRI.
    averageBoundaryError = averageBoundaryError + totalBoundaryError / length(imageLabelCell);
    averageRI = averageRI + sumRI / length(imageLabelCell);
    averageVOI = averageVOI + sumVOI / length(imageLabelCell);
    averageGCE = averageGCE + sumGCE / length(imageLabelCell);

    disp(['Current err:  Boundary  RI  VOI  GCE:']);
    disp([num2str(averageBoundaryError/imageCount) '  ' num2str(averageRI/imageCount) ...
         '  ' num2str(averageVOI/imageCount) '  ' num2str(averageGCE/imageCount)]);
end
A=[BE;RI;VOI;GCE];
st1=strcat(filename,'_I');
xlswrite(st1,A');

averageBoundaryError = averageBoundaryError / imageCount;
averageRI = averageRI / imageCount;
averageVOI = averageVOI / imageCount;
averageGCE = averageGCE / imageCount;


B=[averageBoundaryError;averageRI;averageVOI;averageGCE];
st2=strcat(filename,'_avg_I');
xlswrite(st2,B);
end