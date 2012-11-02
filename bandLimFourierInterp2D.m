%BANDLIMFOURIERINTERP2D 2D Band-Limited Fourier Interpolation
%   fout = bandLimFourierInterp2D(x,f,xout) takes periodic signal f sampled
%   at equispaced grid (x,y) and interpolates to arbitrary set of points 
%   (xout,yout). [x,y] should be a tensor-product grid generated with
%   meshgrid. [xout,yout] can either be a tensor-product grid or a list
%   of points, where xout and yout are each 1D arrays.
%
%   Uses band-limited interpolation formula in tensor product form to avoid 
%   'for' loops. cf. Trefethen, "Spectral Methods in MATLAB", p3.m.
%
%   Author: Derek Steinmoeller, University of Waterloo, 2012.
function fout =  bandLimFourierInterp2D(x,y,f,xout,yout)
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
    
    %make output points into column vectors for now.
    xout=xout(:);
    yout=yout(:);
       
    %build tensor-product interpolation grid & data
    Nout=length(xout);
    %Ndata=length(x);
    [xx,yy,xxout] = meshgrid(x1d,y1d,xout);
    [xx,yy,yyout] = meshgrid(x1d,y1d,yout);
    ff = repmat(f,[1 1 Nout]);
    
    %evaluate all the shifted sinc functions
    fout = ff.*sin(pi*(xxout-xx)/dx)./(pi*(xxout-xx)/dx)...
              .*sin(pi*(yyout-yy)/dy)./(pi*(yyout-yy)/dy);
    
    %double-contract (sum over) them all to get interpolated 
    %data at (xout,yout).
    fout = squeeze(sum(sum(fout,2),1));
    
    %reshape data if we're supposed to
    if TENSOR_PROD_OUTPUT == true
        fout = reshape(fout,szout(1),szout(2));
    end
end