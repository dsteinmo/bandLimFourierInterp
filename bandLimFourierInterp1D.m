%BANDLIMFOURIERINTERP1D 1D Band-Limited Fourier Interpolation
%   fout = bandLimFourierInterp1D(x,f,xout) takes periodic signal f sampled
%   at equispaced grid x and interpolates to arbitrary set of points xout. 
%
%   Uses band-limited interpolation formula in tensor product form to avoid 
%   'for' loops. cf. Trefethen, "Spectral Methods in MATLAB", p3.m.
%
%   Author: Derek Steinmoeller, University of Waterloo, 2012.
function fout =  bandLimFourierInterp1D(x,f,xout)
    %make sure input is column vectors
    x=x(:);
    f=f(:); 
    xout=xout(:);
    %get grid-spacing
    h = abs(x(2)-x(1));
    
    %build tensor-product grid & data
    Nout=length(xout);
    [xx,xxout] = meshgrid(x,xout);
    ff = repmat(f,1,Nout).';
    
    %evaluate all the shifted sinc functions
    ffout = ff.*sin(pi*(xxout-xx)/h)./(pi*(xxout-xx)/h);
    
    %contract (sum over) them all to get interpolated data at xout.
    fout = sum(ffout,2); %sum over columns
end