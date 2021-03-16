
%%
% PA & JN
% Matlab script for
% Stommel/Lorenz basin images with shift
%
% 14th March 2021
%
%
%% Set seed
tic;
clear variables;
rng(1812);

printfigs=true;

%% params for system
p.rho=28;
p.mu=10;
p.beta=8/3;
p.xi1=3;
p.eta1=1;
p.zeta=0.3;
p.A=0.0157;

%% time scaling
p.tau=20;

%% noise
p.sigma=0;

%% shift
shift=-10.0;
% p.rate varies

%% start/end
Tm=-10;
Tp=10;

%% timestep for plotting
h=0.05;
nt=(Tp-Tm)/h+1;

%% start figure count
fignum=0;

%% grid for r
nr=51;

%% number of samples in ensemble
ne=200;

%% Test these values of p.A
As=[0.0135,0.0157];

for i=1:2
    p.A=As(i);
    
    rs=linspace(0.5,1.5,nr);
    nu=zeros(nt,nr);
    
    %% start near Lorenz attractor
    yinit=[1.6,-2.6,26.6,1.7,1.0];
    
    %% simple scan by r
    for ir=1:nr
        p.rate=rs(ir);
        p.shift=shift;
        fprintf('simple scan r=%f\n',p.rate);
        %%
        
        % set initial condition
        yi=yinit+1e-2*rand(size(yinit));
        ti=Tm;
        [tt,y]=runODE(ti,yi,Tp,h,p);
        yy=y(end,:)';
        %% compute indicative norm at each step
        % Psi= "T-S"
        Psi(:,ir)=(y(:,4)-y(:,5)+Lambda(tt,p));
        
    end
    
    fignum=fignum+1;
    %% figs of alpha, epsilon and Gamma
    f2=figure(fignum);
    f2.PaperUnits='centimeters';
    f2.PaperSize=[12 6];
    f2.Units='centimeters';
    f2.InnerPosition=[1 1 12 6];
    clf;
    
    %% plot the time series
    subplot(2,1,1)
    hold on;
    for ir=1:nr
        plot(tt,Psi(:,ir));
    end
    xlabel('$t$','Interpreter','latex');
    ylabel('$\Psi$','Interpreter','latex');
    title(sprintf('$a=%f$',p.A),'Interpreter','latex');
    xlim([tt(1) tt(end)])
    %ylim([0.45 0.55]);
    subplot(2,1,2)
    imagesc(tt,rs,Psi');
    axis xy
    colorbar;
    xlabel('$t$','Interpreter','latex');
    ylabel('$r$','Interpreter','latex');
    
    drawnow();
 hold off;
    
    
    fname=sprintf('SL-BASIN-R-SIMPLE2-%s',date);
    savefigure(fname,fignum,printfigs);
end

toc;

tp=zeros(nr,1);

%% ensemble scan by r
p.A=0.0135;
for ir=1:nr
    p.rate=rs(ir);
    p.shift=shift;
    fprintf('ensemble scan r=%f\n',p.rate);
    %%
    
    % set initial condition
    esc=0;
    parfor i=1:ne
        fprintf('.');
        yi=yinit+1e-2*rand(size(yinit));
        ti=Tm;
        [tt,y]=runODE(ti,yi,Tp,h,p);
        yy=y(end,:)';
        % Psi= "T-S"
        temp=(yy(4)-yy(5)+Lambda(Tp,p));
        % if escaped to off state
        if temp<0.3
            esc=esc+1;
        end
    end
    % probability of tipping
    nesc(ir)=esc;
end
    

fignum=fignum+1;
%% figs of alpha, epsilon and Gamma
f3=figure(fignum);
f3.PaperUnits='centimeters';
f3.PaperSize=[12 6];
f3.Units='centimeters';
f3.InnerPosition=[1 1 12 6];


clf;
rt=rs';
pesc=nesc/ne;
pt=pesc';
rconf = [rt;rt(end:-1:1)];

%% Modified Wald/ Agresti-Coull 95% confidence interval
ptm=(nesc'+2)/(ne+4);
aci=1.96*sqrt(ptm.*(1-ptm)/ne);
pconf = [ptm+aci;ptm(end:-1:1)-aci(end:-1:1)];

p1 = fill(rconf,pconf,'red');
p1.FaceColor = [1 0.8 0.8];      
p1.EdgeColor = 'none';           
hold on
plot(rt,pt,'x-');
hold off
xlabel('$r$','Interpreter','latex');
ylabel('$p_2$','Interpreter','latex');
title(sprintf('$a=%f$',p.A),'Interpreter','latex');
xlim([rt(1) rt(end)]);
ylim([0 1]);

fname=sprintf('SL-SHIFTSCAN-%s',date);
savefigure(fname,fignum,printfigs);

toc;

keyboard;

function dy=FF(yy,t,p)
%% RHS of ODE
% Lorenz + Stommel with shift
%

% shift by Lambda (1,0,0,0,1)
%ll=Lambda(t,p);
temp=(1+tanh(p.rate*t));
ll=p.shift*temp/2;
x1=yy(1)-ll;
x2=yy(2);
x3=yy(3);
x4=yy(4);
x5=yy(5)-ll;

% Lorenz (x1,x2,x3) forcing Stommel (x4,x5)
ydot(1) = p.mu*(x2-x1);
ydot(2) = x1*p.rho-x1*x3-x2;
ydot(3) = x1*x2-p.beta*x3;
ydot(4) = (p.xi1+p.A*x1)-x4*(1+abs(x4-x5));
ydot(5) = (p.eta1+p.A*x1)-x5*(p.zeta+abs(x4-x5));


dy=p.tau*ydot';

end

function [tt,yy]=runODE(ti,yi,tend,h,p)
tspan=ti:h:tend;
[tt,yy]=ode45(@(t,y) FF(y,t,p),tspan,yi);
end

function ll=Lambda(t,p)
temp=(1+tanh(p.rate*t));
ll=p.shift*temp/2;
end

function savefigure(name,number,printfigs)
%
ffname = sprintf('./figs/%s-fig%i', name, number);

h=figure(number);
set(h,'Units','centimeters');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

if(printfigs==true)
    %savefig(ffname);
    print(sprintf('%s.pdf',ffname),'-dpdf','-r0');
end

end