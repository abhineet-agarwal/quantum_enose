clear all
% close all
%%%%%%%%%%%%%%%%%%%%%%%%The role of inelastic scattering in resonant
%%%%%%%%%%%%%%%%%%%%%%%%tunnelling heterostructures%%%%%%%%%%%%%%%%%
%Constants (all MKS, except energy which is in eV)
hbar=1.06e-34;q=1.6e-19;m= 0.067*9.1e-31;IE=(q*q)/(2*pi*hbar);
Ef=0.02;kT=.025;
t0=5.2;
%inputs
a=sqrt((hbar^2)/(2*m*(t0)*q));
% NS=15;NC=30;ND=15;Np=NS+NC+ND;
%Hamiltonian matrix
%NS=15;NC=20;ND=15;Np=NS+NC+ND;UB=0*ones(Np,1);%no barrier
%NS=23;NC=4;ND=23;Np=NS+NC+ND;
%UB=[zeros(NS,1);0.4*ones(NC,1);zeros(ND,1);];%tunneling barrier





%%%%%%%%%%%%%%%%%%%%% one peak in TM vs E %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NS=1;NC=40;ND=1;Np=NS+NC+ND;
Nb=15; Vb=0.6;
UB=[zeros(NS,1);Vb*ones(Nb,1);zeros(NC-2*Nb,1);Vb*ones(Nb,1);zeros(ND,1)];%RT barrier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% original %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NS=15;NC=16;ND=15;Np=NS+NC+ND;
% UB=[zeros(NS,1);0.4*ones(4,1);zeros(NC-8,1);0.4*ones(4,1);zeros(ND,1)];%RT barrier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tf=4.5;


Nf=6;

T=(2*t0*diag(ones(1,Np)))-(t0*diag(ones(1,Np-1),1))-(t0*diag(ones(1,Np-1),-1));
T(Nf,Nf)=2*tf;
T(Nf,Nf-1)=-tf;
T(Nf,Nf+1)=-tf;
% T(Nf-1,Nf)=2*tf;
% T(Nf-1,Nf-1)=-tf;
% T(Nf-1,Nf+1)=-tf;
% T(Nf+1,Nf)=2*tf;
% T(Nf+1,Nf-1)=-tf;
% T(Nf+1,Nf+1)=-tf;
% 
T=T+diag(UB);
%Bias
NV=76;VV=linspace(0,0.3,NV);dV=VV(2)-VV(1);
dE=0.005;
D0=0;
Dnu=[0.1];
Nph=size(Dnu,2);
hnu=[35];%Multiply by dE for actual hnu
Nhnu=1./((exp(dE*hnu./kT))-1);

Dnup=[0.9];
Nphp=size(Dnup,2);
hnup=[18];%Multiply by dE for actual hnu
Nhnup=1./((exp(dE*hnup./kT))-1);

for iV=1:NV
V=VV(iV);mu1=Ef+(V/2);mu2=Ef-(V/2);
U1=V*[.5*ones(1,NS) linspace(0.5,-0.5,NC) -.5*ones(1,ND)];
U1=U1';%Applied potential profile
%Energy grid for Green�s function method
E=[-0.2:dE:0.8];NE=size(E,2);
zplus=i*1e-12;
iV
f1=1./(1+exp((E-mu1)./kT));
f2=1./(1+exp((E-mu2)./kT));


% n=zeros(1,NE);p=zeros(1,NE);
sigin1=zeros(Np,Np,NE);sigout1=zeros(Np,Np,NE);
sigin2=zeros(Np,Np,NE);sigout2=zeros(Np,Np,NE);
sigin=0*ones(Np,Np,NE);sigout=0*ones(Np,Np,NE);
n=zeros(Np,Np,NE);
p=zeros(Np,Np,NE);
gamp=zeros(Np,Np,NE);
gam1=zeros(Np,Np,NE);
gam2=zeros(Np,Np,NE);
G=zeros(Np,Np,NE);

%For infinite 2-D cross-section
%f1=(2*m*kT*q/(2*pi*hbar�2)).*log(1+exp((mu1-E)./kT));
%f2=(2*m*kT*q/(2*pi*hbar�2)).*log(1+exp((mu2-E)./kT));
%Transmission


change=1;
it=1/2; %%%%%for faster convergence
while change>5e-5
for k=1:NE
sig1=zeros(Np);sig2=zeros(Np);sig3=zeros(Np);
ck=1-((E(k)+zplus-U1(1)-UB(1))/(2*t0));ka=acos(ck);
sig1(1,1)=-t0*exp(i*ka);gam1(:,:,k)=i*(sig1-sig1');
ck=1-((E(k)+zplus-U1(Np)-UB(Np))/(2*t0));ka=acos(ck);
sig2(Np,Np)=-t0*exp(i*ka);gam2(:,:,k)=i*(sig2-sig2');
% sig3(Np/2,Np/2)=0;gam3=i*(sig3-sig3');%B�uttiker probe
sigin1(:,:,k)=f1(k)*gam1(:,:,k);sigin2(:,:,k)=f2(k)*gam2(:,:,k);
sigout1(:,:,k)=(1-f1(k))*gam1(:,:,k);sigout2(:,:,k)=(1-f2(k))*gam2(:,:,k);
gamp(:,:,k)=sigin(:,:,k)+sigout(:,:,k);
G(:,:,k)=inv(sparse(((E(k)+zplus)*eye(Np))-T-diag(U1)-sig1-sig2+(i*0.5*gamp(:,:,k))));
A(:,:,k)=(i*(G(:,:,k)-G(:,:,k)'));
n(:,:,k)=(real(G(:,:,k)*((f1(k)*gam1(:,:,k))+(f2(k)*gam2(:,:,k))+sigin(:,:,k))*G(:,:,k)'));
p(:,:,k)=(A(:,:,k)-n(:,:,k));
end
siginnew=D0*n;sigoutnew=D0*p;


for iph=1:Nphp
    ne=zeros(Np,Np,NE);
na=zeros(Np,Np,NE);
pe=zeros(Np,Np,NE);
pa=zeros(Np,Np,NE);
inu=hnup(iph);
if inu<NE
ne(Nf,Nf,[1:NE-inu])=n(Nf,Nf,[inu+1:NE]);%ne=[ne zeros(1,inu)];
na(Nf,Nf,[inu+1:NE])=n(Nf,Nf,[1:NE-inu]);%na=[zeros(1,inu) na];
pe(Nf,Nf,[1:NE-inu])=p(Nf,Nf,[inu+1:NE]);%pe=[pe zeros(1,inu)];
pa(Nf,Nf,[inu+1:NE])=p(Nf,Nf,[1:NE-inu]);%pa=[zeros(1,inu) pa];
siginnew=siginnew+((Nhnup(iph)+1)*Dnup(iph)*ne)+(Nhnup(iph)*Dnup(iph)*na);
sigoutnew=sigoutnew+(Nhnup(iph)*Dnup(iph)*pe)+((Nhnup(iph)+1)*Dnup(iph)*pa);
end
end

for iph=1:Nph
    ne=zeros(Np,Np,NE);
na=zeros(Np,Np,NE);
pe=zeros(Np,Np,NE);
pa=zeros(Np,Np,NE);
inu=hnu(iph);
if inu<NE
    for nn=1:Nf-1
ne(nn,nn,[1:NE-inu])=n(nn,nn,[inu+1:NE]);%ne=[ne zeros(1,inu)];
na(nn,nn,[inu+1:NE])=n(nn,nn,[1:NE-inu]);%na=[zeros(1,inu) na];
pe(nn,nn,[1:NE-inu])=p(nn,nn,[inu+1:NE]);%pe=[pe zeros(1,inu)];
pa(nn,nn,[inu+1:NE])=p(nn,nn,[1:NE-inu]);%pa=[zeros(1,inu) pa];
    end
    for nn=Nf+1:Np
ne(nn,nn,[1:NE-inu])=n(nn,nn,[inu+1:NE]);%ne=[ne zeros(1,inu)];
na(nn,nn,[inu+1:NE])=n(nn,nn,[1:NE-inu]);%na=[zeros(1,inu) na];
pe(nn,nn,[1:NE-inu])=p(nn,nn,[inu+1:NE]);%pe=[pe zeros(1,inu)];
pa(nn,nn,[inu+1:NE])=p(nn,nn,[1:NE-inu]);%pa=[zeros(1,inu) pa];
    end
 
siginnew=siginnew+((Nhnu(iph)+1)*Dnu(iph)*ne)+(Nhnu(iph)*Dnu(iph)*na);
sigoutnew=sigoutnew+(Nhnu(iph)*Dnu(iph)*pe)+((Nhnu(iph)+1)*Dnu(iph)*pa);
end
end


change=sum(sum(sum(abs(siginnew-sigin))));
change=change+sum(sum(sum(abs(sigoutnew-sigout))));
sigin=((1-it)*sigin)+(it*siginnew);
sigout=((1-it)*sigout)+(it*sigoutnew);
end


I1=0;I3=0;Ico=0;Inco=0;I2=0;%Current
for k=1:NE
  I1=I1+real(trace((sigout2(:,:,k)*n(:,:,k))-(sigin2(:,:,k)*p(:,:,k))));%I1=sum(I1);
  I3=I3+real(trace((sigout(:,:,k)*n(:,:,k))-(sigin(:,:,k)*p(:,:,k))));
  I2=I2+real(trace((sigout1(:,:,k)*n(:,:,k))-(sigin1(:,:,k)*p(:,:,k))));
  Ico=Ico+real(trace((sigin2(:,:,k)*G(:,:,k)*gam1(:,:,k)*G(:,:,k)')-(gam2(:,:,k)*G(:,:,k)*sigin1(:,:,k)*G(:,:,k)')));
  Inco=Inco+real(trace((sigin2(:,:,k)*G(:,:,k)*gamp(:,:,k)*G(:,:,k)')-(gam2(:,:,k)*G(:,:,k)*sigin(:,:,k)*G(:,:,k)')));
% I2=real((sigout2.*n)-(sigin2.*p));I2=sum(I2);
% I3=real((sigout.*n)-(sigin.*p));I3=sum(I3);
% I123=dE*[sum(I1) sum(I2) sum(I3)],%Normalized Conductance
%IE=(dE/(V*V))*[sum(E.*I1) sum(E.*I2) sum(E.*I3)],%Normalized Power
%kirchoff=[sum(I123) sum(IE)],%checking for conservation of current and energy current
end
II(iV)=sum(I1)*dE*IE;
II3(iV)=sum(I3)*dE*IE;
II4(iV)=sum(I2)*dE*IE;
IIco(iV)=sum(Ico)*dE*IE;
IInonco(iV)=sum(Inco)*dE*IE;
save(['rtd_flat_foreign_' num2str(hnu*dE) 'eV_Dnu_' num2str(Dnu) 'Nf_' num2str(Nf) 'Dnup_' num2str(Dnup) '_hnuf_' num2str(hnup*dE) 'eV0.mat'])
%save('onlyforeign.mat')
%save(['rtd_flat_foreign_Dnu_' num2str(Dnu) 'Nf_' num2str(Nf) 'Dnup_' num2str(Dnup) '_hnuf_' num2str(hnup*dE) 'eV.mat'])

end

G1=diff(II)./dV;VG=VV([2:NV]);
IETS=diff(G1)./dV;VETS=VV([3:NV]);

%save("-v7","rtd_foreign_flat.mat")
% II2 = smooth(VV,II,0.1,'loess');


[b,a]=butter(2,.1);
II2=filtfilt(b,a,II);
G2=diff(II2)./dV;VG=VV([2:NV]);
g2=filtfilt(b,a,G2);
IETS2=diff(g2)./dV;VETS=VV([3:NV]);
%% IETS2=diff(G2)./dV;VETS=VV([3:NV]);
save(['rtd_flat_foreign_' num2str(hnu*dE) 'eV_Dnu_' num2str(Dnu) 'Nf_' num2str(Nf) 'Dnup_' num2str(Dnup) '_hnuf_' num2str(hnup*dE) 'eV0.mat'])


figure(1)
% % h=plot(VETS,IETS./G1(2:end));
grid on
hold on
%h=plot(VETS,IETS2./g2(2:end),'r');
plot(VETS,diff(G2)./G2(2:end))
ylabel('$(d^2I/dV^2)/(dI/dV)$ $[V^{-1}]$')
xlabel('Bias Voltage [V]')
% saveas(figure(1),['IETS_rtd_foreign_NOfSites' num2str(Np) '.png'])
% saveas(figure(1),['IETS_rtd_foreign_NOfSites' num2str(Np) '.fig'])
% saveas(figure(1),['IETS_rtd_foreign_NOfSites' num2str(Np) '.eps'], 'epsc2')
% % saveas(figure(1),['IETS_rtd_NOfSites' num2str(Np) '.png'])
% % saveas(figure(1),['IETS_rtd_NOfSites' num2str(Np) '.fig'])
% % saveas(figure(1),['IETS_rtd_NOfSites' num2str(Np) '.eps'], 'epsc2')

% % saveas(figure(1),['IETS_rtd_sat_NOfSites' num2str(Np) '.png'])
% % saveas(figure(1),['IETS_rtd_sat_NOfSites' num2str(Np) '.fig'])
% % saveas(figure(1),['IETS_rtd_sat_NOfSites' num2str(Np) '.eps'], 'epsc2')




figure(2)
plot(VV,II,'r')
hold on
grid on
% plot(VV,II2)
ylabel('Current [A]')
xlabel('Bias Voltage [V]')
% saveas(figure(2),['IV_rtd_foreign_NOfSites' num2str(Np) '.png'])
% saveas(figure(2),['IV_rtd_foreign_NOfSites' num2str(Np) '.fig'])
%% saveas(figure(2),['IV_rtd_foreign_NOfSites' num2str(Np) '.eps'], 'epsc2')
%% % saveas(figure(2),['IV_rtd_NOfSites' num2str(Np) '.png'])
%% % saveas(figure(2),['IV_rtd_NOfSites' num2str(Np) '.fig'])
%% % saveas(figure(2),['IV_rtd_NOfSites' num2str(Np) '.eps'], 'epsc2')
%
%% saveas(figure(2),['IV_rtd_sat_NOfSites' num2str(Np) '.png'])
%% saveas(figure(2),['IV_rtd_sat_NOfSites' num2str(Np) '.fig'])
%% saveas(figure(2),['IV_rtd_sat_NOfSites' num2str(Np) '.eps'], 'epsc2')



% % % figure(3)
% % % % h=plot(VETS,IETS);
% % % hold on
% % % h=plot(VETS,IETS2);
% % % grid on
% % % % legend(Legend)
% % % saveas(figure(3),['IETS_NN_rtd_np_NOfSites' num2str(Np) '.png'])
% % % saveas(figure(3),['IETS_NN_rtd_np_NOfSites' num2str(Np) '.fig'])
% % % saveas(figure(3),['IETS_NN_rtd_np_NOfSites' num2str(Np) '.eps'], 'epsc2')
% % saveas(figure(3),['IETS_NN_rtd_NOfSites' num2str(Np) '.png'])
% % saveas(figure(3),['IETS_NN_rtd_NOfSites' num2str(Np) '.fig'])
% % saveas(figure(3),['IETS_NN_rtd_NOfSites' num2str(Np) '.eps'], 'epsc2')

% % saveas(figure(3),['IETS_NN_sat_np_NOfSites' num2str(Np) '.png'])
% % saveas(figure(3),['IETS_NN_sat_np_NOfSites' num2str(Np) '.fig'])
% % saveas(figure(3),['IETS_NN_sat_np_NOfSites' num2str(Np) '.eps'], 'epsc2')

% figure(4)
% hold on
% plot(VV,II3,VV,II4)
% legend('Phonon contact','second contact')
% % plot(VV,II)
% hold on
% grid on
% % plot(VV,II2)
% ylabel('Current [A]')
% xlabel('Bias Voltage [V]')
% % % saveas(figure(2),['IV_rtd_np_NOfSites' num2str(Np) '.png'])
% % % saveas(figure(2),['IV_rtd_np_NOfSites' num2str(Np) '.fig'])
% % % saveas(figure(2),['IV_rtd_np_NOfSites' num2str(Np) '.eps'], 'epsc2')
% saveas(figure(4),['I3V_rtd_NOfSites' num2str(Np) '.png'])
% saveas(figure(4),['I3V_rtd_NOfSites' num2str(Np) '.fig'])
% % saveas(figure(2),['IV_rtd_NOfSites' num2str(Np) '.eps'], 'epsc2')
% 
% 
% figure(5)
% hold on
% plot(VV,-IIco,VV,-IInonco,VV,II2)
% legend('Coherent','Non-coherent','Total')
% % plot(VV,IIincoh)
% hold on
% grid on
% % plot(VV,II2)
% ylabel('Current [A]')
% xlabel('Bias Voltage [V]')
% % % saveas(figure(2),['IV_rtd_np_NOfSites' num2str(Np) '.png'])
% % % saveas(figure(2),['IV_rtd_np_NOfSites' num2str(Np) '.fig'])
% % % saveas(figure(2),['IV_rtd_np_NOfSites' num2str(Np) '.eps'], 'epsc2')
% saveas(figure(5),['IVall_rtd_NOfSites' num2str(Np) '.png'])
% saveas(figure(5),['IVall_rtd_NOfSites' num2str(Np) '.fig'])
% % saveas(figure(2),['IV_rtd_NOfSites' num2str(Np) '.eps'], 'epsc2')