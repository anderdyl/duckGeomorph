
clear
cd('/home/dylananderson/projects/duckGeomorph/')

load('/media/dylananderson/Elements1/NC_climate/Nags_Head_SLPs_2degree_memory.mat','slp_mem','time','m_regional','n_regional','XRsq','YRsq','X_in','Y_in')
Xs = min(min(XRsq)):2:max(max(XRsq));
Ys = min(min(YRsq)):2:max(max(YRsq));
X_B = XRsq(:);
Y_B = YRsq(:);

[XR,YR] = meshgrid(Xs,Ys);

for qq = 1:length(X_in)
    sea_nodes(qq) = find(XR == X_in(qq) & YR == Y_in(qq));
end

time = datenum(time(:,1),time(:,2),time(:,3));
data = slp_mem;
data(:,71) = data(:,70);
% a = 0.1:.1:10;
% for i = 0:199
%   data((i+1),:) = sin(pi*(0.5*a + 0.1*i)) + 5*(rand(1,100)-0.5);
% end
% data = (data - ones(200,1)*mean(data));
% time = [0:199];

[lamda, loadings, pcs, per] = ceof(data,10);

figure
pcolor(data)
shading interp

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

mode = 4;
figure
subplot(2,2,1)
plot(S(:,mode).*sqrt(lamda(mode)),'.')
ylabel('Spatial amplitude')
subplot(2,2,2)
plot(theta2(:,mode),'.')
ylabel('Spatial phase')
subplot(2,2,3)
plot(time,Rt(:,mode)./sqrt(lamda(mode)),'o')
datetick
ylabel('Temporal amplitude')
subplot(2,2,4)
plot(time,phit2(:,mode),'o')
datetick('x')
ylabel('Temporal phase')


cd('/media/dylananderson/Elements1/shusin6_contents/codes/')
figure
subplot(2,2,[1 3])
[m_regional, n_regional] = size(XR);
load('paleta2.mat');
load('mycmap_col.mat')
wt = ones(m_regional*n_regional,1)*NaN;
temp = S(:,mode).*sqrt(lamda(mode));
wt(sea_nodes) = temp;
reshaped = reshape(wt,m_regional,n_regional);
    iso1_N = min(temp):(max(temp)-min(temp))/20:max(temp);
    m_proj('miller','lat',[-55 60],'lon',[-90 30]);
    m_coast('patch',[0.5 0.5 0.5]);
    m_grid('box','fancy','xtick',[],'ytick',[])
    hold on
    [cs,h] = m_contourf(XR-360,YR,reshaped,iso1_N,'linewidth',0.5,'linecolor','none');

subplot(2,2,[2 4])
wt2 = ones(m_regional*n_regional,1)*NaN;
temp2 = theta2(:,mode);
wt2(sea_nodes) = temp2;
reshaped = reshape(wt2,m_regional,n_regional);
    iso1_N = min(temp2):(max(temp2)-min(temp2))/20:max(temp2);
    m_proj('miller','lat',[-55 60],'lon',[-90 30]);
    m_coast('patch',[0.5 0.5 0.5]);
    m_grid('box','fancy','xtick',[],'ytick',[])
    hold on
    [cs,h] = m_contourf(XR-360,YR,reshaped,iso1_N,'linewidth',0.5,'linecolor','none');

figure
mode = 10
wtreal = ones(m_regional*n_regional,1)*NaN;
tempreal = eofreal(:,mode);
wtreal(sea_nodes) = tempreal;
wtimag = ones(m_regional*n_regional,1)*NaN;
tempimag = eofimag(:,mode);
wtimag(sea_nodes) = tempimag;

reshapedReal = reshape(wtreal,m_regional,n_regional);
reshapedImag = reshape(wtimag,m_regional,n_regional);

iso1_N = min(temp2):(max(temp2)-min(temp2))/20:max(temp2);

    m_proj('miller','lat',[-55 60],'lon',[-90 30]);
    m_coast('patch',[0.5 0.5 0.5]);
    m_grid('box','fancy','xtick',[],'ytick',[])
    hold on
    
    %[cs,h] = m_contourf(XR-360,YR,reshaped,iso1_N,'linewidth',0.5,'linecolor','none');
    m_quiver(XR-360,YR,reshapedReal,reshapedImag)

