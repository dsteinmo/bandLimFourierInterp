%BANDLIMFOURIERINTERP1D 1D Band-Limited Fourier Interpolation
%   fout = bandLimFourierInterp1D(x,f,xout) takes periodic signal f sampled
%   at equispaced grid x and interpolates to arbitrary set of points xout. 
%
%   Uses band-limited interpolation formula in tensor product form to avoid 
%   'for' loops. cf. Trefethen, "Spectral Methods in MATLAB", p3.m.
%
%   Author: Derek Steinmoeller, University of Waterloo. Copyright (C) 2012.
%   Thanks: Mike Dunphy (improvements to speed/memory efficiency)
%           Chris Subich (help with correcting periodicity bugs)
%
%   This program is free software; you can redistribute it and/or
%   modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation; either version 2
%   of the License, or (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA

function fout =  bandLimFourierInterp1D(x,f,xout)

    %make sure input is column vectors
    x=x(:);
    f=f(:); 
    xout=xout(:);
    
    %get grid-spacing, number of nodes, and domain length.
    dx = abs(x(2)-x(1));
    Nx = length(x);
    L = dx*length(x);
    
    %build tensor-product grid & data
    Nout=length(xout);
    [xx,xxout] = meshgrid(x,xout);
    ff = repmat(f,1,Nout).';
    
    %depdening on the parity of the original signals number of nodes
    %choose the appropriate periodic sinc.
    if mod(Nx,2)==1
        psinc = @(x) bsinc(pi*x,dx,L)*(dx/L);
    else
        psinc = @(x) bsinc(pi*x,dx,L).*cos(pi*x/L).*(dx/L);
    end
    
    %compute kernel function and calculate the contribution
    %due to each grid point
    kernel = psinc(xxout-xx);
    ffout = ff.*kernel;
    
    %contract (sum over) grid points to get interpolated data at xout.
    fout = sum(ffout,2); %sum over columns
end
% Helper function: Baseline sinc -- returns sin(x/dx)/sin(x/L), since
% the special case x == 0 gives a naive implementation NaN issues
function z = bsinc(x,dx,L)
    z = zeros(size(x));
    z(x==0) = L/dx;
    z(x~=0) = sin(x(x~=0)/dx)./sin(x(x~=0)/L);
    return
end