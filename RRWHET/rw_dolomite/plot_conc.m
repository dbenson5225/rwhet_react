clear all
%geom=zeros(5,6);
%concfile=fopen('dolomite.con');
vel=2.4e-5;Diff=1.2e-7;x_0=0.05;
concfile='dolomite.con';
geom=dlmread(concfile);
nspec=geom(1,7)
times=(size(geom,1)-7)/(2+2*nspec);
mom=zeros(times,4);
nx=geom(1,1);ny=geom(1,2);nz=geom(1,3);
dx=geom(1,4);dy=geom(1,5);dz=geom(1,6);
dirs=0;
for i=1:3
    if geom(1,i)>1.5
        dirs=dirs+1;
    end
end
nxy=nx*ny;
maxc=max(geom(11,:))
%domain=squeeze(zeros(nx,ny,nz));
%domain=zeros(nx,ny,nz,nspec);

for l=1:times
    domain=zeros(nx,ny,nz,nspec);
    adv=(2+2*nspec)*(l-1);
    t=geom(adv+8,1);
    nboxes=geom(adv+8,2);
% convert box numbers to i,j,k
    for n=1:nboxes
        ijk=geom(adv+9,n);
        kbox=floor((ijk-1)/nxy)+1;
        jbox=floor((ijk-1-(kbox-1)*nxy)/nx)+1;  
        ibox=ijk-(kbox-1)*nxy-(jbox-1)*nx;
        for m=1:nspec;
          domain(ibox,jbox,kbox,m)=geom(adv+9+2*m,n);  % puts conc in right box
        end
    end
% Plot in the right number of dimensions
  if dirs==1
%    cpos=domain(domain>0);
    for m=1:nspec
      xplot=dx*(1:nx)';
      cplot=domain(:,:,:,m);
      xplot=xplot(cplot>0); cplot=cplot(cplot>0)
      figure(m);plot(xplot,cplot);     
      axis([0 nx*dx 0 maxc])
      mean=sum(xplot.*cplot)/sum(cplot);
      var=sum(xplot.*xplot.*cplot)/sum(cplot)-mean.^2;
      mom(l,m,1)=t; mom(l,m,2)=mean; mom(l,m,3)=var; mom(l,m,4)=2*Diff*t;
      mom(l,m,5)=sum(cplot)*dx;  
    end
      pause
  end
  if dirs==2
    contourf(squeeze(domain));  % not sure which dimension is singleton
  end
  if dirs==3  
  end
end   % time loop

figure; plot(mom(:,1,1),mom(:,1,2),'-o'); hold on
plot(mom(:,1,1), x_0+vel*mom(:,1,1));
figure; plot(mom(:,1,1),mom(:,1,3),'-o'); hold on
plot(mom(:,1,1),mom(:,2,3),'-o')
plot(mom(:,1,1),mom(:,1,4));
figure(); plot(mom(:,1,1),mom(:,1,5));
title(['Time = ',num2str(t),', Mean = ',num2str(mean),', Var = ',num2str(var)])
%    outfile=['plots/plot',num2str(t)]
%    print('-f11',outfile,'-dpdf');
%    pause
    