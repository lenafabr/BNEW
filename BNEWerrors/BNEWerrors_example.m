% example script for calculating errors in parameters from BNEW analysis

% ---------
%% Single persistent random walk drift with tau=30
% ---------

% load in data matrices generated by the fortran code
% results are saved in a .mat file for quicker loading later
datamats1 = covartxt2mat('covar.t30');

D = 3; % diffusion coefficient
locE = 1; %localization error
nvals = 3:15; % wavelet spans to use
ntrack = 100;

kscl=0.74; % maximal k value to go to as a factor of n

% calculate errors as a function of gamma
ng=50;
gamlist = linspace(0,3,ng);
bias = zeros(ng,3); errvals = zeros(ng,3); toterrvals = zeros(ng,3);

for gc = 1:ng
    gam = gamlist(gc);
    
    [bias(gc,:),errvals(gc,:),toterrvals(gc,:),trscl,MSDrscl] = getFitParamErrors(nvals,gam,D,locE,ntrack,datamats1,'kscl',kscl);
end

%% plot results
plot(gamlist,toterrvals(:,1)/D,'b',gamlist,toterrvals(:,2),'k',gamlist,toterrvals(:,3),'r','LineWidth',1.5)
hold all
plot(gamlist,bias(:,1)/D,'b--',gamlist,bias(:,2),'k--',gamlist,bias(:,3),'r--','LineWidth',1.5)
hold off
xlabel('gamma')
ylabel('bias and error')
legend('alpha', 'D','locE')

% ---------
%% Two independent persistent random walk drifts with tau=30, tau=10
% ---------
% load in data matrices generated by the fortran code
% results are saved in a .mat file for quicker loading later
datamats1 = covartxt2mat('covar.t30');
datamats2 = covartxt2mat('covar.t10');
couplemat = dlmread('couple.t30.t10.txt');

% calculate errors as a function of D
gam = [1,1]
ng=50;
Dlist = linspace(0.1,5,ng);
bias = zeros(ng,3); errvals = zeros(ng,3); toterrvals = zeros(ng,3);
kscl=0.74;

for gc = 1:ng
    D = Dlist(gc);
    
    [bias(gc,:),errvals(gc,:),toterrvals(gc,:),trscl,MSDrscl] = getFitParamErrors_2exp(nvals,gam,D,locE,ntrack,datamats1,datamats2,couplemat,'kscl',kscl);
end

%% plot results
plot(Dlist,toterrvals(:,1)./Dlist','b',Dlist,toterrvals(:,2),'k',Dlist,toterrvals(:,3),'r','LineWidth',1.5)
hold all
plot(Dlist,bias(:,1)./Dlist','b--',Dlist,bias(:,2),'k--',Dlist,bias(:,3),'r--','LineWidth',1.5)
hold off
xlabel('D')
ylabel('bias and error')
legend('alpha', 'D','locE')

ylim([-0.2,0.4])