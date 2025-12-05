function data = gauss_smooth(data,FWHM,pad)
    %function data = gauss_smooth(data,FWHM)
    %smooth in xy plane with a Gaussian Kernel
    %replaces input image with smoothed image
    
    %data = input volume time series 
    %xdim = x dimension
    %ydim = y dimension
    %zdim = z dimension
    %tdim = t dimension
    %data_out = smoothed image
    %
    % data = gauss_smooth(data,2);
    
    
    if (~exist('pad','var')),
        pad = 0;
    end;
    
    if (pad > 0),
        pad1 = zeros(pad,1);
        data = [pad1;flatrow(data);pad1];
    end;
    
    tdim = length(data);
    %create Gaussian Kernel
    s  = [FWHM];
    
    s  = s/sqrt(4*log(2));				% FWHM -> Gaussian parameter
    
    x  = round(6*s(1));
    kdimx = 2*x+1;
    x = [-x:x];
    x  = exp(-(x).^2/(2*(s(1)).^2)); 
    x  = x/sum(x);
    
    gauss = zeros(kdimx,1);
    for i = 1:kdimx,
      gauss(i) = sqrt(x(i));
    end;
    gauss = gauss/max(gauss);
    gauss = gauss/sum(gauss);
    [val ind] = max(gauss);
    %ind = round(ind/2);
    data = mean(data)+conv(data-mean(data),gauss);
    %data = data(ind-1:tdim+ind-1);
    data = data(ind:tdim+ind-1);
    
    if (pad > 0)
        data = data(pad+1:(end-pad));
    end;
end

