function [lamda, loadings, pcs, per] = ceof(data, nkp);

%  a = 0.1:.1:10;
%  for i = 0:199
%    data((i+1),:) = sin(pi*(0.5*a + 0.1*i)) + 5*(rand(1,100)-0.5);
%  end
%  data = (data - ones(200,1)*mean(data));
  
  if nargin < 2; nkp = 10; end;
  if nargin < 1; error('Need to input data'); end;

  [ntim, npt] = size(data);

  disp('Calculating hilbert matrix')
  %data = data + j * hilbert(data);
  data2 = hilbert(data);
  disp('Done with hilbert matrix, calculating covariance matrix')
  c = data2' * data2 / ntim;
  c2 = data2' * data2;
  disp('Covariance matrix computed, starting eig(c)')

  [loadings, lamda] = eig(c);
  l = diag(lamda);

  [lamda,k] = sort(l'); loadings = loadings(:,k);
  lamda      = fliplr(lamda);
  loadings   = fliplr(loadings);
  loadings = loadings(:,1:nkp);
  per = real(lamda*100/sum(lamda));
  pcs = data2 * loadings;

%  loadings = loadings(:,1:100);
%  pcs = data * loadings(:,1:10);