function out = normMeanFrame(mf,excludeFreq)
    if nargin < 2 || isempty(excludeFreq)
        excludeFreq = 15;
    end

    I3=fft2(mf);
    I3=fftshift(I3);
    I4=log(1+abs(I3));
    [x y] = meshgrid([1:length(I3(1,:))]-length(I3(1,:))./2+0.5, ...
        [1:length(I3(:,1))]-length(I3(:,1))./2+0.5);
    d2c = sqrt(x.^2 + y.^2);
    mask = d2c<excludeFreq;
    I5=I3.*mask;
    I4=log(1+abs(I5));
    I6=ifft2(ifftshift(I5));
    lpf = abs(I6);
    
    out = mf-lpf;
end