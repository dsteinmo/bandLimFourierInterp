%unit test for bandLimtFourierInterp2D - modified from
% Trefethen, "Spectral Methods in MATLAB", p3.m.
clear
close all;

%make sure that the functions we're calling are visible to matlab
addpath ../

%set error tolerance
errtol = 1e-3;

%build grid
h = .5; 
xmax = 10; 
xmin = -xmax;
x = -xmax:h:xmax;

%pick some interpolation points
xout = randn(1,200)*xmax;
xout = xout(xout <= xmax & xout >= xmin); %make sure we're _interpolating_
%xout = sort(xout,'ascend');

%set function to interpolate, and do interpolation to above points
v = sech(x).^2;
p = bandLimFourierInterp1D(x,v,xout);

%check error - should get better if we refine the grid
err = norm(sech(xout(:)).^2-p(:),2)/length(p);
if err < errtol
    disp(['Test of bandLimFourierInterp1D.m PASSED with err=' num2str(err)]);
else
    disp(['Test of bandLimFourierInterp1D.m FAILED with err=' num2str(err)]);
end

%%plot results if you want
% figure(1); clf;
% plot(x,v,'.-',xout,p,'.r'); 
% legend('original function','interpolated values');
% axis([xmin xmax -.1 1.1]); 

% %test 2 - an aperiodic function -> this is expected to fail due to gibbs
% oscillations, since input data must be periodic
% v = tanh(x);
% xout=x+h/2;
% p = bandLimFourierInterp1D(x,v,xout);
% 
% err = norm(tanh(xout(:))-p(:),2)/length(p);
% if err < errtol
%     disp(['Test of bandLimFourierInterp1D.m (aperiodic) PASSED (somehow?!) with err=' num2str(err)]);
% else
%     disp(['Test of bandLimFourierInterp1D.m (aperiodic) FAILED (as expected) with err=' num2str(err)]);
% end
% 
% figure(2); clf;
% plot(x,v,'.-',xout,p,'.r'); 
% legend('original function','interpolated values');
% axis([xmin xmax -1.1 1.1]); 