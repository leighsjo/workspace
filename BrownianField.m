function outField = BrownianField(H,imSize,outputFlag, isPlotted)
%% simulate Fractional Brownian field on unit disk, with Hurst parameter 'H';
%  Note that the covariance function is isotropic, see reference below.
% INPUTS:
%        - 'H' is the Hurst parameter of the Gaussian process
%        - 'imSize' is the number of grid points in row and column: 
%            Both 'm' and 'n' should be a power of 2;
%            if not, they are set to be by e.g. n=2^ceil(log2(n)); 
%        - 'outputFlag':
%           1) 0 (default): get absolute field
%           2) 1          : get real part of the field
%           3) 2          : get imaginary part of the field
%        - 'isPlotted'
%           1) 0 (default)
%           2) 1: plot the field
% OUTPUT:
%          - two statistically independent fields 'field1' and 'field2'
%            over unit disk; if not output requested, then function
%            outputs a figure of one of the fields
%          - vectors 'tx' and 'ty' so that the field is plotted via
%            surf(tx,ty,field1,'EdgeColor','none')
%
% Example:
%  [field1,field2,tx,ty]=Brownian_field(.9,2^10);
%   surf(tx,ty,field2,'EdgeColor','none'),colormap bone
%% Reference:
% Kroese, D. P., & Botev, Z. I. (2015). Spatial Process Simulation.
% In Stochastic Geometry, Spatial Statistics and Random Fields(pp. 369-404)
% Springer International Publishing, DOI: 10.1007/978-3-319-10064-7_12
if (H>1)||(H<0) % Hurst parameter error check
    error('Hurst parameter must be between 0 and 1')
end
if nargin < 4
    isPlotted = 0;
    if nargin < 3
        outputFlag = 0;
        if nargin < 2
            imSize = 2^8; % default value of points
        else
            imSize =2.^ceil(log2(imSize));
        end
    end
end
if length(imSize) == 1
    [m,n] = deal(imSize);
elseif length(imSize) == 2
    [m,n] = deal(imSize(1),imSize(2));
end
R=2; % [0,R]^2 grid, may have to extract only [0,R/2]^2
tx=(1:n)/n*R; ty=(1:m)/m*R; % create grid for field
Rows=zeros(m,n);
for i=1:n
    for j=1:m % rows of blocks of cov matrix
        Rows(j,i)=rho([tx(i),ty(j)],[tx(1),ty(1)],R,2*H);
    end
end
BlkCirc_row=[Rows, Rows(:,end-1:-1:2);
    Rows(end-1:-1:2,:), Rows(end-1:-1:2,end-1:-1:2)];
% compute eigen-values
lam=real(fft2(BlkCirc_row))/(4*(m-1)*(n-1));
lam=sqrt(lam);
% generate field with covariance given by block circular matrix
Z=complex(randn(2*(m-1),2*(n-1)),randn(2*(m-1),2*(n-1)));
F=fft2(lam.*Z);
F=F(1:m,1:n); % extract sub-block with desired covariance
[~,~,c2]=rho([0,0],[0,0],R,2*H);
switch outputFlag
    case 0 
        outField = getOutputField(F,tx,ty,c2,@abs);
    case 1
        outField = getOutputField(F,tx,ty,c2,@real);
    case 2
        outField = getOutputField(F,tx,ty,c2,@imag);
end

if isPlotted
    surf(tx,ty,outField,'EdgeColor','none')
%     title('Fractional Gaussian Field on Unit Disk')
    colormap bone
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out,c0,c2]=rho(x,y,R,alpha)
% embedding of covariance function on a [0,R]^2 grid
if alpha<=1.5 % alpha=2*H, where H is the Hurst parameter
    beta=0;c2=alpha/2;c0=1-alpha/2;
else % parameters ensure piecewise function twice differentiable
    beta=alpha*(2-alpha)/(3*R*(R^2-1)); c2=(alpha-beta*(R-1)^2*(R+2))/2;
    c0=beta*(R-1)^3+1-c2;
end
% create continuous isotropic function
r=sqrt((x(1)-y(1))^2+(x(2)-y(2))^2);
if r<=1
    out=c0-r^alpha+c2*r^2;
elseif r<=R
    out=beta*(R-r)^3/r;
else
    out=0;
end
end
function outField = getOutputField(F,tx,ty,c2,fH)
outField=fH(F);
outField=outField-outField(1,1); % set field zero at origin
% make correction for embedding with a term c2*r^2
outField=outField + kron(ty'*randn,tx*randn)*sqrt(2*c2);
% [X,Y]=meshgrid(tx,ty);
% outField((X.^2+Y.^2)>1)=nan;
% outField=outField(1:n/2,1:m/2);
% varargout{1} = outField;
% if nargout >1
%     tx=tx(1:n/2);ty=ty(1:m/2);
%     varargout{2} = tx;
%     varargout{3} = ty;
% end
end