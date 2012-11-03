%BANDLIMFOURIERINTERP2D 2D Band-Limited Fourier Interpolation
%   fout = bandLimFourierInterp2D(x,f,xout) takes periodic signal f sampled
%   at equispaced grid (x,y) and interpolates to arbitrary set of points 
%   (xout,yout). [x,y] should be a tensor-product grid generated with
%   meshgrid. [xout,yout] can either be a tensor-product grid or a list
%   of points where xout and yout are each 1D arrays.
%
%   fout = bandLimFourierInterp2D(x,f,xout,maxMem) is the same as above but
%   with optional argument maxMem that allows the specification of maximum
%   amount of memory (in bytes) allowed to be allocated by temporary 3D
%   tensor-product arrays. Default is 200 MB.
%
%   Uses band-limited interpolation formula in tensor product form to avoid 
%   'for' loops. cf. Trefethen, "Spectral Methods in MATLAB", p3.m.
%
%   Author: Derek Steinmoeller, University of Waterloo, 2012.
function fout =  bandLimFourierInterp2D(x,y,f,xout,yout,maxMem)
 
    global xx yy xxout yyout ff

    if ~exist('maxMem','var')
        maxMem = 209715200;%default = 200 MB
    end

    %get 1D versions of x & y-grids
    x1d = x(1,:);
    y1d = y(:,1);
    Nx = length(x1d);
    Ny = length(y1d);

    %get grid-spacing
    dx = abs(x1d(2)-x1d(1));
    dy = abs(y1d(2)-y1d(1));
    
    %determine if we need to return a 1D or 2D array.
    szout = size(xout);
    if szout(1) > 1 && szout(2)> 1
        TENSOR_PROD_OUTPUT = true;
    else
        TENSOR_PROD_OUTPUT = false;
    end
    
    %make output points into column vectors for rest of calculations
    xout=xout(:);
    yout=yout(:);
    Nout=length(xout);
    
    %calculate acceptable size of 3D arrays
    sizeOf3D =maxMem/5;               %5 3D arrays
    maxEntries3D = sizeOf3D/8;        %populated by 8-byte entries.
    %so third dimension should be at most this big:
    zlen = floor(maxEntries3D/Nx/Ny);
    
    if zlen == 0
        error('input data too big or maxMem too small. Try making maxMem larger (default is 209715200 bytes).');
    end
    
    %if we can fit the entire calculation in memory 
    %at the same time, then proceed as normal (as in older version)
    if Nout < zlen
        %allocate memory for 3D arrays
        ff = zeros(Nx,Ny,Nout);
        
        %build tensor-product interpolation grid & data
        [xx,yy,xxout] = meshgrid(x1d,y1d,xout);
        [xx,yy,yyout] = meshgrid(x1d,y1d,yout);
        
        fout = doSincInterp(f,dx,dy);        
    else
        %if we can't fit it into memory, split it up into blocks of length
        %zlen and process them sequentially with a for-loop.
        Nblocks = ceil(Nout/zlen);
        fout = zeros(Nout,1);
        
        for n=1:Nblocks
            %build tensor-product interpolation grid & data
            if n==Nblocks
                %last block could be smaller
                range = (zlen*(n-1)+1):Nout;
            else
                range = (zlen*(n-1)+1):zlen*n;
            end
            [xx,yy,xxout] = meshgrid(x1d,y1d,xout(range));
            [xx,yy,yyout] = meshgrid(x1d,y1d,yout(range));
            
            fout(range) = doSincInterp(f,dx,dy);    
        end
    end
    
    %reshape data if we're supposed to
    if TENSOR_PROD_OUTPUT == true
        fout = reshape(fout,szout(1),szout(2));
    end
end

%Helper function that actually does the interpolation - to be called 
%repeatedly with pre-allocated global 3D arrays.
function fout = doSincInterp(f,dx,dy)
    global xx yy xxout yyout ff
    
    szout = size(xx);
    Nout = szout(3);
    %evaluate all the shifted sinc functions
    
    ff = repmat(f,[1 1 Nout]);
    fout = ff.*sin(pi*(xxout-xx)/dx)./(pi*(xxout-xx)/dx)...
              .*sin(pi*(yyout-yy)/dy)./(pi*(yyout-yy)/dy);
    
    %double-contract (sum over) them all to get interpolated 
    %data at (xout,yout).
    fout = squeeze(sum(sum(fout,2),1));
end