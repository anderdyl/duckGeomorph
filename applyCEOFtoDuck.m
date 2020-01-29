clear
data = ncread('/home/dylananderson/projects/alllinesSouth.nc','alllines');
time = ncread('/home/dylananderson/projects/alllinesSouth.nc','t');

time = double(time)+datenum(1981,10,1);
demeandata = data'-mean(data');

[lamda, loadings, pcs, per] = ceof(demeandata, 200);

%loadings = loadings.*sqrt(lamda');
%pcs = pcs./sqrt(lamda);
% EOFs
eofreal = real(loadings(:,1:10));
eofimag = imag(loadings(:,1:10));
% PCs
pcreal = real(pcs(:,1:10));
pcimag = imag(pcs(:,1:10));
% Spatial amplitude
S = power(loadings.*conj(loadings),0.5);
% Spatial phase
theta = atan2(eofimag,eofreal);
theta2 = theta.*180./pi;
% Temporal amplitude
Rt=power(pcs.*conj(pcs),0.5);
% Temporal  phase
phit = atan2(pcimag,pcreal);
phit2 = phit.*180./pi;

mode = 1;
figure
subplot(2,2,1)
plot(S(:,mode),'.')
ylabel('Spatial amplitude')
subplot(2,2,2)
plot(theta2(:,mode),'.')
ylabel('Spatial phase')
subplot(2,2,3)
plot(time,Rt(:,mode),'o')
datetick
ylabel('Temporal amplitude')
subplot(2,2,4)
plot(time,phit2(:,mode),'o')
datetick('x')
ylabel('Temporal phase')



%   Outputs:
%   lamda    = eigenvalues (should be pretty close to real)
%   loadings = First 10 Complex Loadings (row dimension)
%   pcs      = First 10 Complex Principal Components (column dim)
%   per      = percent variance explained (real)
% 
%   Inputs:
%   data     = data, so that size(data) = [ntime nspace]
%   nkp      = number of modes to output (default = 10);
% 
%   Note:  pcs can be found by performing the following:
%     pcs = data * loadings(:,1:10);
% 
%   Normalization is such that:
%     loadings' * loadings = diag(ones(1,npt));
%     pcs' * pcs = diag(lamda);
% 
%   For display purposes, the following patterns and time
%     series go together:
%          real(loadings) goes with real(pcs);
%          imag(loadings) goes with imag(pcs) = hilbert(real(pcs))
% 
%   Also, one can divide the pcs by sqrt(lamda) and multiply the
%     loadings by sqrt(lamda) to get actual amplitudes.  Recall,
%     std(real(pcs)) should equal std(imag(pcs)).

