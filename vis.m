%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% visualization of the instabile mode %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all
addpath('freefem_matlab_octave_plot-public/release-v2.0/demos/ffmatlib')
addpath('Data')
addpath('Eigenvalues')
addpath('Eigenvectors')
addpath('Flows')

%loading data from ff++
%load mesh
[p,b,t,nv,nbe,nt,labels]=ffreadmesh('sudden_exp.msh');
%load connectivity
VPh=ffreaddata('sudden_exp_VPh.txt');
%load baseflow
[Ub] = ffreaddata('Vectored_baseflow.txt');
[VhP2,ub]=ffvectorget({'P2','P2','P1'}, VPh, Ub, 1);
[~,vb]=ffvectorget({'P2','P2','P1'}, VPh, Ub, 2);
[VhP1,pb]=ffvectorget({'P2','P2','P1'}, VPh, Ub, 3);
%load eigevalues
id = fopen('eig.txt');
d = textscan(id,'%f %f');
e = cell2mat(d);
fclose('all');

%baseflow visualization
figure(1)
ffpdeplot(p,b,t,'VhSeq',VhP2,'XYData',ub,'Mesh','off','MColor','b','Boundary','off','ColorMap',jet);
title('Baseflow x component')
xlabel('x')
ylabel('y')
%number of eigenvalues calculated
nEV = length(e(:,1));

%load eigenvectors and build eigenvalues in cartesian form
la = zeros(nEV,1);
up = zeros(length(ub),nEV);
vp = zeros(length(vb),nEV);
pp = zeros(length(pb),nEV);

for i =1:nEV
    la(i) = e(i,1)+1i*e(i,2);
    name = sprintf('vectorized_eigenvectors_%i.txt',i);
    [Up] = ffreaddata(name);
    [~, up(:,i)]=ffvectorget({'P2','P2','P1'}, VPh, Up, 1);
    [~, vp(:,i)]=ffvectorget({'P2','P2','P1'}, VPh, Up, 2);
    [~, pp(:,i)]=ffvectorget({'P2','P2','P1'}, VPh, Up, 3);
end

%complex plane plot
figure(2)
hold on
set(gca, 'XAxisLocation', 'origin')
set(gca, 'YAxisLocation', 'origin')
plot(la, 'o');
plot(conj(la), 'o');
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
grid on
xlim([min(real(la))*0.9 max(real(la))*1.1])
ylim([-max(imag(la))*1.1 max(imag(la))*1.1])

%selecting the eigenvalue withe the greatest real part and the relative
%modes

[~,i] = max(real(la));

lai = la(i);
upi = up(:,i);
vpi = vp(:,i);
ppi = pp(:,i);

%linearized flow
A = 0.25; %sufficently large amplitude to visualize the perturbation

temp = 0:0.01:100;
ind=2669;
y=zeros(1,length(temp));
for i = 1: length(temp)
    v = vb+A*real(vpi*exp(lai*temp(i)));
    y(i)=v(ind);
end
figure(3)
plot(temp,y);
xlabel('t')
ylabel('v')
grid on

T = 50;
figure(10)

utemp=real(up);
vtemp=real(vp);
ffpdeplot(p,b,t,'VhSeq',VhP2,'XYData',utemp,'Mesh','off','MColor','b','Boundary','off','ColorMap',jet);


