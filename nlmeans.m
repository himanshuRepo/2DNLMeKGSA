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
% File purpose: Computes the non-local means of an image.
%===================================================================================

function [sig_est]=nlmeans(sig_orig,t,f)

[M,N]=size(sig_orig);

%===================================================================================
% parameters: warning - change at your own risk - especially modifying "W",
%      "wnd_size" or "V" may significantly affect computation time - only change
%      these parameters when you know what you are doing.
% 2nd note - changing some of the parameters (e.g. wnd_size, weighting_function_name)
%      may also require to use a different parameter value for "h")
%===================================================================================

sigma = 0.2;               % noise standard deviation for the simulation
correlated_noise = 0;     % 0 - white Gaussian noise
                          % 1 - correlated Gaussian noise                                                   

W = f;                    % half window size (actual window size = 2*W+1) 
                          % so W=4 corresponds to 9x9 blocks
if correlated_noise==0,
h = 2.08*sigma;           % "bandwidth" parameter
else
h = 2.48*sigma;           % "bandwidth" parameter    
end;

wnd_size = t;             % half size of the _search_ window
                          % note: normally I use wnd_size=15 (31x31 window)
                          % or larger, but not here, because matlab is
                          % simply too slow:)
weighting_function_name = 'MODIFIED_BISQUARE'; % the weighting function being used
vector_based_filter = 1;  % use the vector based NLMeans filter
klt_postprocessing = 1;   % use KLT-based postprocessing to denoise patches
                          % with a low number of similar candidates
V = 2;                    % half size of the vector neighborhood
                          % (we use V<W mainly to limit memory usage)


if correlated_noise == 0,                          
    sig = sig_orig + sigma*randn(M,N);  % test signal with white noise
else    
    [img_noise,PS]=psdnoise(M,N,13.2,1);
    sig = sig_orig + sigma*img_noise;  % test signal with correlated noise
end;
sig=sig_orig;
                          
%===================================================================================
% IMPORTANT REMARK: with klt_postprocessing=1, if the bandwidth parameter is chosen
% too small, there will be a lot of noise remaining in the NLMeans-denoised image.
% the KLT-denoising filter will compensate for that and will remove the remainder of
% the noise, however the results then correspond to a purely local filter. The 
% (standalone) KLT-filter has a comparable performance as wavelet-based denoising 
% filters such as, e.g. BLS-GSM.

% Next - pad the input signal (to deal with boundary effects - i.e. when the
% patches at the borders of the signal are being compared). 
sig = bound_extension(sig,W,W,'mirror');
N=N+2*W;
M=M+2*W;
Bs=(2*W+1)^2; % block size

%% PART I - initialization & preprocessing
switch weighting_function_name
    case 'DEFAULT' % LECLERC
        weighting_function = @(r) exp(-r.^2/(h*h)); 
    case 'TUKEY'
        weighting_function = @(r) (r<=h) .* (1-r.*r)/(h*h);
    case 'HUBER'
        weighting_function = @(r) 1./max(1,r/h);
    case 'ANDREWS'
        weighting_function = @(r) (r==0) .* 1 + (r~=0 && r<=h) .* sin(pi*r/h)./(r/h);
    case 'CAUCHY'
        weighting_function = @(r) 1/(1+r.*r/(h*h));
    case 'GEMANMCCLURE'
        weighting_function = @(r) 1./((1+r.*r/(h*h))*(1+r.*r/(h*h)));
    case 'LECLERC'
        weighting_function = @(r) exp(-r.^2/(h*h)); 
    case 'FAIR'
        weighting_function = @(r) 1./(1+r/h);
    case 'BISQUARE'
        weighting_function = @(r) (r<h).*(1-(r/h).^2);
    case 'MODIFIED_BISQUARE'
        weighting_function = @(r) (r<h).*(1-(r/h).^2).^8;
    case 'LOGISTIC'
        weighting_function = @(r) (r==0) .* 1 + (r~=0) .* tanh(r/h)/(r/h);
    case 'TALWAR'
        weighting_function = @(r) (r<h).*1;     
    case 'BLUE'
        weighting_function = @(r) (r<h)./(h*h) + (r>=h)./(r.*r);
end;

% the KLT-postprocessing filter is designed to work in combination with the
% vector-based NLMeans filters
if klt_postprocessing && ~vector_based_filter
    error('The KLT-based postprocessing can only be applied if vector_based_filter==1!');
end;

% we're dealing with correlated noise => compute the prewhitened signal
if correlated_noise == 1,        
    sig_prewhit=real(ifft2(fft2(sig(W+1:end-W,W+1:end-W))./max(1e-4,abs(PS)).^0.5));
    sig_prewhit = bound_extension(sig_prewhit,W,W,'mirror');

else
    sig_prewhit=sig; % no prewhitening necessary
end;

%% PART II - NLMeans filtering
fprintf('NLMeans filtering');

if vector_based_filter == 0,
    % the final image starts with all pixel intensities initialized to 0.
    sig_est = zeros(size(sig));
else
    % we're estimating (2V+1 x 2V+1) local neighborhood vectors. Therefore
    % we need a bit more memory... :(
    sig_est = zeros(size(sig,1),size(sig,2),(2*V+1)^2);
end;

% two accumulation buffers (one for the weights, a second for the local noise variance)
weights_accum = zeros(size(sig));
sq_weights_accum = zeros(size(sig));      % for variance computation


for md = -wnd_size:wnd_size,
    for nd = -wnd_size:wnd_size,          % the shift (Delta_i in our LNLA paper)        
        if md>0 || (md==0 && nd>0),
            fprintf('.');
            
            % compute the feasible ranges for m and n (to stay inside the signal boundaries)
            m_min = max(min(W-md,M-W),W+1);
            m_max = min(max(M-W-1-md, W),M-W-1);
            n_min = max(min(W-nd, N-W),W+1);
            n_max = min(max(N-W-1-nd, W),N-W-1);
            
            if n_min>n_max || m_min>m_max,
                continue;
            end;             
            
            % here follows the moving average filter implementation from Section 5.2 from
            % our LNLA paper... The purpose here is to show how the algorithm works (as in the
            % C++ version). I would like to remark that the for-loops are of course not optimal 
            % for MATLAB. For a MATLAB implementation, you may want to use either separable 
            % filtering or recursive filtering functions from MATLAB instead (to do so, set
            % use_moving_average = 0).
            use_2dmoving_average = 0;
            
            if use_2dmoving_average,           
                % initialization
                square_diff = zeros(1,N);
                square_diff2 = zeros(1,N);
                ssd = zeros(M,N);
                ssd_cum = 0;
                m = m_min;
                m1 = m+md;            
                
                % moving average warming up-phase I
                for n=(n_min-W):(n_min+W-1)
                    square_diff(1+n) = sum(sum((sig_prewhit(1+m+(-W:W),1+n)-sig_prewhit(1+m1+(-W:W),1+n+nd)).^2));
                    ssd_cum = ssd_cum + square_diff(1+n);
                end;    

                % moving average warming up-phase II
                for n=n_min:n_max,
                    n1 = n + nd;
                    square_diff(1+n+W) = sum(sum((sig_prewhit(1+m+(-W:W),1+n+W)-sig_prewhit(1+m1+(-W:W),1+n1+W)).^2));
                    ssd_cum = ssd_cum + square_diff(1+n+W);
                    ssd(1+m,1+n) = ssd_cum;
                    ssd_cum = ssd_cum - square_diff(1+n-W);
                end;    

                % running ssd phase
                for m=(m_min+1):m_max
                    m1=m+md;
                    ssd_cum2=0;
                    for n=(n_min-W):(n_max+W)
                        diff=sig_prewhit(1+m-W-1,1+n)-sig_prewhit(1+m1-W-1,1+n+nd);
                        square_diff(1+n)=square_diff(1+n)-diff*diff;
                    end;
                    for n=(n_min-W):(n_min+W-1)
                        % compute the squared difference
                        diff=sig_prewhit(1+m+W,1+n)-sig_prewhit(1+m1+W,1+n+nd);
                        square_diff2(1+n)=diff*diff;
                        square_diff(1+n)=square_diff(1+n)+square_diff2(1+n);
                        ssd_cum2=ssd_cum2+square_diff(1+n);
                    end;
                    for n=n_min:n_max
                        n1=n+nd;
                        % compute the squared difference
                        diff=sig_prewhit(1+m+W,1+n+W)-sig_prewhit(1+m1+W,1+n1+W);
                        square_diff2(1+n+W)=diff*diff;
                        square_diff(1+n+W)=square_diff(1+n+W)+square_diff2(1+n+W);
                        ssd_cum2=ssd_cum2+square_diff(1+n+W);
                        ssd(1+m,1+n)=ssd_cum2;
                        ssd_cum2=ssd_cum2-square_diff(1+n-W);
                    end;
                end;    
                
            else
                
                range_m = 1+(m_min:m_max);
                range_n = 1+(n_min:n_max);
            
                % implementation using separable filtering ...
                ssd = zeros(M,N);
                ssd(range_m,range_n) = conv2 (conv2 ((sig_prewhit(range_m,range_n) - sig_prewhit(range_m+md,range_n+nd)).^2, ...
                        ones(1, 2*W+1), 'same'), ones(2*W+1, 1), 'same');                
            end;


            % compute the weights and use the accumulation technique as described
            % in Section 5.1 of our LNLA paper
            range_m = 1+(m_min:m_max);
            range_n = 1+(n_min:n_max);
            weights = weighting_function(sqrt(ssd(range_m,range_n)/Bs));
            
            if vector_based_filter == 0,
                % accumulate signal + weight
                sig_est(range_m,range_n)            = sig_est(range_m,range_n) + weights .* sig(range_m+md,range_n+nd);
                weights_accum(range_m,range_n)      = weights_accum(range_m,range_n) + weights;
                sq_weights_accum(range_m,range_n)   = sq_weights_accum(range_m,range_n) + weights.^2;

                % exploiting the weight symmetry
                sig_est(range_m+md,range_n+nd)       = sig_est(range_m+md,range_n+nd) + weights .* sig(range_m,range_n);
                weights_accum(range_m+md,range_n+nd) = weights_accum(range_m+md,range_n+nd) + weights;
                sq_weights_accum(range_m+md,range_n+nd) = sq_weights_accum(range_m+md,range_n+nd) + weights.^2;
            else
                % weight accumulation remains the same...
                weights_accum(range_m,range_n)       = weights_accum(range_m,range_n) + weights;
                weights_accum(range_m+md,range_n+nd) = weights_accum(range_m+md,range_n+nd) + weights;
                sq_weights_accum(range_m,range_n)    = sq_weights_accum(range_m,range_n) + weights.^2;
                sq_weights_accum(range_m+md,range_n+nd) = sq_weights_accum(range_m+md,range_n+nd) + weights.^2;

                for m=-V:V,
                    for n=-V:V
                        index=(V+m)*(2*V+1)+V+n+1;
                        sig_shifted = circshift(sig,[m n]);
                        sig_est(range_m,range_n,index)=sig_est(range_m,range_n,index)+weights.*sig_shifted(range_m+md,range_n+nd);
                        sig_est(range_m+md,range_n+nd,index)=sig_est(range_m+md,range_n+nd,index)+weights.*sig_shifted(range_m,range_n);
                    end;
                end;
            end;
        elseif md==0 && nd == 0,
            % experimental feature - changing the weight for the (0,0)-displacement seems to improve
            % the denoising performance
            weight = 0.01*weighting_function(0);
            weights_accum = weights_accum + weight*ones(size(weights_accum));
            sq_weights_accum = sq_weights_accum + weight^2*ones(size(weights_accum));

            if vector_based_filter == 0,
                sig_est = sig_est + weight*sig;
            else
                for m=-V:V,
                    for n=-V:V
                        index=(V+m)*(2*V+1)+V+n+1;
                        sig_est(:,:,index) = sig_est(:,:,index) + weight*circshift(sig,[m n]);
                    end;
                end;
            end;
        end;
    end;
end;

fprintf('\n');

%% PART III - KLT-based postprocessing

if klt_postprocessing,

    % locally adaptive KLT-based denoising
    B = 16; % block size (to reduce memory usage)

    % compute the local noise variance in the image
    local_noise_var = sigma.^2 * sq_weights_accum./max(1e-10,weights_accum).^2;
    
    % a possible refinement for the case of correlated noise would be 
    % to also take the noise PSD into account in the KLT denoising. Here
    % we simply assume that the noise PSD is flat after NLMeans denoising, 
    % which is not necessarily true. Although for most types of correlated
    % noise this does not pose too many problems...

    fprintf('KLT based postprocessing...\n');

    for m=1:B:size(sig_est,1)
        for n=1:B:size(sig_est,2)
            block_rm = m:min(m+B-1,M);
            block_rn = n:min(n+B-1,N);
            sample = sig_est(block_rm, block_rn,:);
            for k=1:size(sample,3) % normalize sample by the weight
                sample(:,:,k) = sample(:,:,k) ./ max(1e-10,weights_accum(block_rm, block_rn));
            end;
            sample = reshape(sample,[size(sample,1)*size(sample,2), size(sample,3)]);
            y_mean = kron(mean(sample,2),ones(1,size(sample,2)));
            Cy = cov(sample - y_mean);
            [U,Sig_Y_diag,dummy] = svd(Cy);
            Sig_Y_diag = diag(Sig_Y_diag).';

            u = (sample - y_mean) * U; % forward KLT

            for m1=block_rm,
                for n1=block_rn,
                    Sig_X_diag = max(1e-2,Sig_Y_diag-local_noise_var(m1,n1));
                    index=(m1-m)*length(block_rn)+ n1-n+1;
                    u(index,:) = u(index,:) .* Sig_X_diag ./ (Sig_X_diag + local_noise_var(m1,n1));
                end;
            end;

            y = u * U' + y_mean; % backward KLT
            y = reshape(y,[length(block_rm),length(block_rn),size(sample,2)]);
            
            for k=1:size(y,3) % unnormalize estimated sample value
                y(:,:,k) = y(:,:,k) .* max(1e-10,weights_accum(block_rm, block_rn));
            end;
            
            % ... and update the estimated signal vectors
            sig_est(block_rm,block_rn,:) = y;
        end;
    end;

end;

%% PART IV - aggregation

if vector_based_filter == 0,
    % normalize by dividing by the weights
    sig_est = sig_est ./ max(1e-20,weights_accum);    
else
    % in the vector-based case, we only accumulated the weights for the central pixel
    % so for. With a simple convolution, we can obtain the block-accumulated weights
    weights_accum = max(1e-20,conv2(weights_accum,ones(2*V+1,2*V+1),'same'));
    sig_final = zeros(size(sig));
    % finally, average over all shifted estimated (this technique is in fact similar
    % to DWT cycle spinning)
    for m=-V:V,
        for n=-V:V
            index=(V+m)*(2*V+1)+V+n+1;
            sig_final = sig_final + circshift(sig_est(:,:,index),-[m n]) ./ weights_accum;
        end;
    end;
    sig_est = sig_final;
    clear sig_final;
end;

% Unpad the output signal (i.e. crop the mirrored borders)
sig_est = sig_est(W+1:end-W,W+1:end-W);
sig = sig(W+1:end-W,W+1:end-W);

%%
% And plot the results
end
% figure,
% subplot(121),imshow(sig,[]),title(sprintf('Noisy image PSNR=%f dB', GetPSNR(sig,sig_orig)));
% subplot(122),imshow(sig_est,[]),title(sprintf('Denoised image PSNR=%f dB', GetPSNR(sig_est,sig_orig)));
