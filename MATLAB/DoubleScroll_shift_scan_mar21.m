
%%
% PA & JN
% Matlab script for
% Double Scroll with rate-dependent shift
%
% March 2021
%
%
%% Set seed
tic;
clear variables;
rng(1811);

printfigs=true;
savedata=false;

% output fig names
fname=sprintf('DSSHIFTSCAN-%s',date);

%% params for functions
p.a=9;
p.b=14;

%% noise
p.eta=0;

%% p.sigma and p.rate vary

%% total duration
total=80;

%% timestep
h=1e-2;
ni=floor(total/h)+1;

%% set up figure numbering
fignum=1;


%% parameter scan
npar=31;
ratemin=0.1;
ratemax=3.1;

%% probability of escape
rates=zeros(npar,1);
nesc=zeros(npar,1);

%% Threshold for escape
ythresh=6;

%% trial repeats
nreps=400;

%% initial value
yi=[0,0,0];

yf=zeros(nreps,3);

% consider two values of shift
for sh=1:2
    if sh==1
        ratemin=0.1;
        ratemax=10.1;
        p.sigma=1;
    else
        ratemin=0.1;
        ratemax=3.1;
        p.sigma=2;
    end
    for i=1:npar
        p.rate=ratemin+(ratemax-ratemin)*(i-1)/(npar-1);
        rates(i)=p.rate;
        f1=figure(fignum);
        %f1.PaperUnits='centimeters';
        %f1.PaperSize=[12 10];
        f1.Units='centimeters';
        f1.InnerPosition=[1 1 14 10];
        clf;
        endnorm=zeros(nreps,1);
        for j=1:nreps
            k=1;
            t=-total/2;
            y=zeros(ni,3);
            tt=zeros(ni,1);
            rns=randn(size(tt));
            yy=yi+[2*(rand-0.5),0.5*(rand-0.5),0.5*(rand-0.5)];
            y(k,:)=yy;
            tt(k)=t;
            while (t<total/2)
                % Euler-Maruyama step
                fftemp=FF(yy,t,p);
                ynew=yy+h*fftemp+...
                    p.eta*sqrt(h)*rns(k)*yy;
                t=t+h;
                k=k+1;
                y(k,:)=ynew;
                tt(k)=t;
                yy=ynew;
            end
            yf(j,:)=yy;
            
            
            for ii=1:3
                subplot(3,1,ii)
                plot(tt,y(:,ii));
                hold on;
                xlabel('t');
                ylabel(sprintf('x_%d',ii));
                xlim([tt(1) tt(end)])
                %ylim([0.45 0.55]);
            end
            endnorm(j)=norm([3*(yy(2)-Lambda(t,p)),yy(3)]);
        end
        % number of escapes
        nesc(i)=sum(endnorm>ythresh);
        
        subplot(3,1,1);
        title(sprintf('r=%f p_2=%f sigma=%f',p.rate,nesc(i)/nreps,p.sigma));
        hold off;
        
        drawnow();
        
        if sum(p.rate==[0.1 0.2 1.0 2.9])>0
            savefigure(fname,fignum,printfigs);
            fignum = fignum+1;
        end
        
    end
    
    %savefigure(fname,fignum,printfigs);
    
    
    fignum=fignum+1;
    %% figs of alpha, epsilon and Gamma
    f2=figure(fignum);
    clf;
    %f2.PaperUnits='centimeters';
    %f2.PaperSize=[12 10];
    f2.Units='centimeters';
    f2.InnerPosition=[1 1 14 10];
    
    
    clf;
    rt=rates';
    pesc=nesc/nreps;
    pt=pesc';
    rconf = [rt rt(end:-1:1)];
    %% Modified Wald/ Agresti-Coull interval
    ptm=(nesc'+2)/(nreps+4);
    aci=1.96*sqrt(ptm.*(1-ptm)/nreps);
    pconf = [ptm+aci ptm(end:-1:1)-aci(end:-1:1)];
    
    p1 = fill(rconf,pconf,'red');
    p1.FaceColor = [1 0.8 0.8];
    p1.EdgeColor = 'none';
    
    hold on
    plot(rt,pt,'x-');
    hold off
    title(sprintf('sigma=%f',p.sigma));
    xlabel('r');
    ylabel('p_2');
    xlim([rt(1) rt(end)]);
    ylim([0 1]);
    
    savefigure(fname,fignum,printfigs);
    fignum = fignum+1;
    
end

toc;

keyboard

%savefigure(fname,fignum,printfigs);

function dy=FF(y,t,p)
%% RHS of ODE
% Chua RHS
%

% shift y1
xx=y(1)-Lambda(t,p);
yy=y(2)-Lambda(t,p);
zz=y(3);

dy(1) = p.a*(yy-Phi(xx));
dy(2) = xx-yy+zz;
dy(3) = -p.b*yy;

end

function tt=Phi(xx)
% Hirsch Smale Devaney nonlinearity
%
tt=xx^3/16-xx/6;
end

function ll=Lambda(t,p);
ll=p.sigma*(1+tanh(p.rate*t))/2;
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