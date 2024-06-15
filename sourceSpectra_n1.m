
% -------------------------------------------------------------------------
%   PROGRAM TO MAKE ESTIMATE SOURCE PARAMETERS
%           Dated    : 26th AUGUST 2013 
% -------------------------------------------------------------------------
%       LOADING THE INPUT WAVEFORM DATA IN SAF FORMAT
%--------------------------------------------------------------------------
clc;clear all
for i=1:3 
flg=true;   % set the logic variable to start the loop
while flg
  inputfile=uigetfile('*___001_BH_2N_SAF','MultiSelect','on');
textfilename=inputfile;
fid = fopen(textfilename,'r');
data = textscan(fid,'%f%f%f','delimiter','\t','headerlines', 12);
data{1};

DATA=[data{1} data{2} data{3}];
fclose(fid);
l=length(DATA(:,1));
z=DATA(:,1);%Vertical Component
z=z-mean(z);
y=DATA(:,2);%NS Component
y=y-mean(y);
x=DATA(:,3);%EW Component
x=x-mean(x);
DATA(:,1)=x;DATA(:,2)=y;DATA(:,3)=z;
% -------------------------------------------------------------------------
%   STEP-I          ROTATION Horizontal component(N & E) CLOCKWISE  
% -------------------------------------------------------------------------
azi=122
theta1=deg2rad(azi);
x1=(x.*cos(theta1))-(y.*sin(theta1));
y1=(x.*sin(theta1))+(y.*cos(theta1));
z1=z;
rotated(:,1)=x1;rotated(:,2)=y1;rotated(:,3)=z1;
save rotated.txt rotated* -ascii;
rdata=rotated;
% -------------------------------------------------------------------------
%       INSTRUMENT RESPONSE CORRECTION
%__________________________________________________________________________
%INSTRUMENT CONSTANT= [{SENSITIVITY (V/m/s)}/{GEN CONSTT (V/Counts)}]*100
for i=1:3
A=3.9844e-07; %A=1.6e-07; % Converts counts to cm /sec
CX=rdata(:,i)*A;
p=[-.707+.707j ; -.707-.707j ; -62.3816+135.392j ; -62.3816-135.392j ; -350 ; -75 ];
%p=[-4.44+4.44j;-4.44-4.44j];%CMG-40T1
%p=[-0.036+.038j ; -.036-.038j ; -222+222j ; -222-222j ];%KS2000 General
z=[0;0];
[b,a]=zp2tf(p,z,1.7071e-009);
sys=tf(b,a);
bode(sys);
VDATA(:,i)=filter(b,a,CX);
save VDATA.txt VDATA* -ascii;
end
VX=VDATA(:,1);
%--------------------------------------------------------------------------
plot(VX);title('Select S Arrival');
ylabel('Velocity(cm/sec)');grid on;
% -------------------------------------------------------------------------
PM=ginput(2);
Mp=round(PM(1));
Ms=round(PM(2));
M=Ms;
tsp=(Ms-Mp)/100;%30.12;
tsr=tsp+1.38*tsp;
gr=tsp*8*100000;
SR=100;%input('Enter the sampling rate, SR: ');
SP=1/SR; % Sampling Period
N=1024;%input('Enter the number of samples,N : ');
ln=length(VX);
t=(0:SP:(ln-1)*SP)';
tmin=0;tmax=max(t);
%--------------------------------------------------------------------------
%                       WINDOW FUNCTION 
%--------------------------------------------------------------------------
w=kaiser(N);%     COSINE TAPER  (20%)tukeywin(N,0.2);%
for i=1:3
wdata(:,i)=VDATA(M:M+N-1,i).*w;
end
save wdata.txt wdata* -ascii;
% -------------------------------------------------------------------------
MX=max(abs(VX));
scrsz = get(0,'ScreenSize');% To assign the size of plotting window
figure ('position',[scrsz]);
subplot(3,1,1);plot(t,VX);title('SH Component');
xlabel('time(sec)');ylabel('cm/sec');axis([0 tmax -MX MX])
grid on;hold on

lt=max(MX);ht=2*lt;
rectangle('Position', [M*SP -lt N*SP ht],'FaceColor','g'...
               ,'EdgeColor','red','LineWidth',1);
plot(t(M+1:M+N),VX(M+1:M+N),'r');

% -------------------------------------------------------------------------
saveas(gcf,'plotwdata.fig');
saveas(gcf,'plotwdata.jpg');
%saveas(gcf,'plotwdata.jpg');pause;hold off;close
% -------------------------------------------------------------------------
n=N/2;  
f=(1:n)*(SR/N);
datav(:,1)=f;
datad(:,1)=f;
dataa(:,1)=f;
for i=1:3
    df=wdata(:,i);
    b=fft(df);
    b=abs(b(1:n))./n;
    b=1.33*b;
    %q=5.4.*(1+f/0.3).^3.5./(f/0.3).^2;%82.*f.^1.12;%
    %atnn=exp(-pi*f*tsr./q);%a=exp(nu./de)
    bb1=b'%.*atnn;
    bb=bb1*gr;
    datav(:,i+1)=bb;
    datad(:,i+1)=bb./(2*pi*f);
    dataa(:,i+1)=bb.*(2*pi*f);
end
save datav.txt datav* -ascii;
save datad.txt datad* -ascii;
save dataa.txt dataa* -ascii;
%--------------------------------------------------------------------------
 %            PROGRAM TO FIT SOURCE MODEL 
%  Brune model with high-cut fmax model fitted to source spectrum
%--------------------------------------------------------------------------
%clc;clear all;
for i=1:3
load datav.txt;
fv=datav(:,1);
vel=datav(:,4);
displ=vel./(2*pi*fv).^1;acc=vel.*(2*pi*fv).^1;tecc=vel.*(2*pi*fv).^2.5;
loglog(fv,acc,'m');grid on;title ('SH - Spectrum');
xlabel('Frequency (Hz)'); ylabel('Acceleration Spectrum');grid on;hold on
% =========================================================================
%                ESTMATION OF initial values of fc and fmax
% -------------------------------------------------------------------------
v1=vel;
index1=find(abs(v1)==max(abs(v1)));
fc=datav(index1,1);

t1=tecc;
index2=find(abs(t1)==max(abs(t1)));
fmax=datav(index2,1);
 
avacc=acc(index1:index2);
omg=mean(avacc)/(2*pi*fc).^2;
% -------------------------------------------------------------------------
%                   Estimation of slope above fmax
%--------------------------------------------------------------------------
n=length(fv);
px1=fv(index2);
py1=acc(index2);
px2=fv(n-20);
py2=acc(n-20);
slop=log10(py1/py2)/(log10(px1/px2));
p=abs(slop);
% -------------------------------------------------------------------------
f=fv;lf=length(f);

%       Evaluation of kappa according to least error fit
%--------------------------------------------------------------------------
     kappa=0.040:0.040:0.08;
    lkappa=length(kappa);
    for i=1:lkappa
        kappa1=kappa(i);
    outkappa(1,i)=kappa1;
    spektraf=omg./((1.+(fv./fc).^2).*exp(pi*fv*kappa1));
    outspek(:,i)=spektraf;
    defr(i)=mean(abs(displ-spektraf));
    outkappa(2,i)=defr(i);
    end
    index5=find(abs(defr)==min(abs(defr)));
    kappa=kappa(index5);

%             Estimation of best - fc
% -------------------------------------------------------------------------
    fcw=fv(1:index2);
    lfc=length(fcw);
    for i=1:lfc
        fc1=fcw(i);
    outfc(1,i)=fc1;
    spektraf=omg./((1.+(fv./fc1).^2).*(1.+(fv./fmax).^p));
    outspek(:,i)=spektraf;
    defr(i)=mean(abs(displ-spektraf));
    outfc(2,i)=defr(i);
end
 plot(outfc(1,:),outfc(2,:));
index3=find(abs(defr)==min(abs(defr)));
fc=outfc(1,index3);
% -------------------------------------------------------------------------
%             Estimation of best - fmax
% -------------------------------------------------------------------------
    fmaxw=fv(index1+1:lf);
    lfmax=length(fmaxw);
    for i=1:lfmax
        fmax1=fmaxw(i);
    outfmax(1,i)=fmax1;
    spektrafm=omg./((1.+(fv./fc).^2).*(1.+(fv./fmax1).^p));
    outspekmax(:,i)=spektraf;
    defrm(i)=mean(abs(displ-spektrafm));
    outfmax(2,i)=defrm(i);
end
% plot(outfmax(1,:),outfmax(2,:));
index4=find(abs(defrm)==min(abs(defrm)));
fmax=outfmax(1,index4);

%--------------------------------------------------------------------------
%                Calculation of Source Parameters
%--------------------------------------------------------------------------
m0=(4*pi*2.8*(3.6*100000)^3*omg)/(2*0.6);
mw=(0.6667*log10(m0))-10.7;
r=(2.34*(3.6*100000)/(2*pi*fc))/100;
sd=(7*m0)/(16*(r*100)^3)*0.000001;
m0=m0/10000000;%dyne-cm to Nm
sd=sd/10;%bars to MPa
%SAVE SOURCE PARAMETERS
save ("out.txt","m0","omg","fc","fmax","p","mw","r","sd","kappa","-ascii");hold on;
A =load('out.txt');
for K = 1: 9: numel(A); 

C=transpose(A);
C2=reshape(C,1,9);
%disp(C2(K : min(K:K+9))); 
%C2=unique(C1,'row', 'stable');
save("out_new.txt","C2",'-ascii','-append');
end
%--------------------------------------------------------------------------
f=fv;lf=length(f);
% -------------------------------------------------------------------------
%           Plot fitted fmax spektra
% -------------------------------------------------------------------------
subplot(1,2,1)
loglog(fv,acc,'m');grid on;%title ('SH - Spectrum');
xlabel('frequency (Hz)'); ylabel( 'Acceleration Spectrum');grid on;hold on
spektra=((2*pi*f).^2)*omg./((1.+(f./fc).^2).*(1.+(f./fmax).^p));% Brune Model with high cut
spektra1=((2*pi*fc).^2)*omg./((1.+(fc./fc).^2).*(1.+(fc./fmax).^p));
spektra2=((2*pi*fmax).^2)*omg./((1.+(fmax./fc).^2).*(1.+(fmax./fmax).^p));

%               kappa
spektrak=((2*pi*f).^2)*omg./((1.+(f./fc).^2).*exp(pi*fv*kappa));
spektrak1=((2*pi*fc).^2)*omg./((1.+(fc./fc).^2).*exp(pi*fc*kappa));

                 %plotting for fmax
 loglog(f,spektra,'b','LineWidth',2);grid on;hold on
 plot(fc,spektra1,'ro', 'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','g');
 plot(fmax,spektra2,'ro', 'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','g');
 plot(fc,spektra1,'r+', 'MarkerSize',25); plot(fmax,spektra2,'r+', 'MarkerSize',25);
 %              plotting for kappa
 loglog(f,spektrak,'k--','LineWidth',2);grid on;
 plot(fc,spektrak1,'ro', 'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','g');
 
%plot(fc,spektrak1,'r+', 'MarkerSize',25); 
mostr=num2str(m0, '%3.1e\n');fcstrr=num2str(r, '%3.1f');
omstr=num2str(omg, '%3.1e\n');mstrp=num2str(p,'%3.1f');
mstr1=num2str(mw,'%3.1f');fcstr1=num2str(fc, '%3.1f');
mstr2=num2str(sd,'%3.1f');fcstr2=num2str(fmax, '%3.1f');
strk=num2str(kappa, '%4.3f');

	
tx=fv(1);
% sl=length(spektra);
ty=min(spektra);
text(tx,ty*3,[ ' M_0= ',mostr,' Nm;','  {\it f}_c= ',fcstr1,' Hz;',...
    ' {\it f}_m_a_x= ', fcstr2,' Hz;',' p = ' ,mstrp,  ' sd= ' ,mstr2,]);	
 
%text(tx,ty*3,[' M_0= ',mostr,' Nm;','  {\it f}_c= ',fcstr1,' Hz;',...
%   ' {\it f}_m_a_x= ', fcstr2,' Hz;',' p = ' ,mstrp,],'HorizontalAlignment','left',..
 %  'BackgroundColor',[1 1 .9],'EdgeColor','r','LineWidth',1,'FontSize',12);
   
text(tx,ty*3,['M_0= ',mostr,' Nm;','  {\it f}_c= ',fcstr1,' Hz;',...
    ' {\it f}_m_a_x= ', fcstr2,' Hz;',' p = ' ,mstrp,' k= ' ,strk,],'HorizontalAlignment','left',...
    'BackgroundColor',[1 1 .9],'EdgeColor','r','LineWidth',1,'FontSize',12); 
text(tx,ty,[' M_w= ',mstr1,';','{\it r} = ',fcstrr,' m;', '\Delta\sigma=',mstr2,' MPa'],...
    'HorizontalAlignment','left','BackgroundColor',[1 1 .9],'EdgeColor','r',...
        'LineWidth',1,'FontSize',12);
		
    saveas(gcf,'fmaxaccautofit.fig');saveas(gcf,'acc_spectra.jpg');

%	saveas(gcf,'fmaxaccautofit.fig');saveas(gcf,'acc_spectra.jpg');pause;close
% -------------------------------------------------------------------------
%                     PLOTTING OF DISPLACEMENT SPECTRA
% -------------------------------------------------------------------------
subplot(1,2,2)
loglog(fv,displ,'m');grid on;%title ('SH - Spectrum');
xlabel('frequency (Hz)'); ylabel( 'Displacement Spectrum');grid on;hold on;

% -------------------------------------------------------------------------
spektrd=((2*pi*f).^0)*omg./((1.+(f./fc).^2).*(1.+(f./fmax).^p));% Brune Model with high cut
spektrd1=((2*pi*fc).^0)*omg./((1.+(fc./fc).^2).*(1.+(fc./fmax).^p));
spektrd2=((2*pi*fmax).^0)*omg./((1.+(fmax./fc).^2).*(1.+(fmax./fmax).^p));
%               kappa
spektrak=((2*pi*f).^0)*omg./((1.+(f./fc).^2).*exp(pi*fv*kappa));
spektrak1=((2*pi*fc).^0)*omg./((1.+(fc./fc).^2).*exp(pi*fc*kappa));

% plot of fmax
 loglog(f,spektrd,'b','LineWidth',2);grid on;
 plot(fc,spektrd1,'ro', 'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','g');
 plot(fmax,spektrd2,'ro', 'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','g');
 plot(fc,spektrd1,'r+', 'MarkerSize',25); plot(fmax,spektrd2,'r+', 'MarkerSize',25);
 
 %              plotting for kappa
 loglog(f,spektrak,'k--','LineWidth',2);grid on;
 plot(fc,spektrak1,'ro', 'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','g');
 plot(fc,spektrak1,'r+', 'MarkerSize',25); 
 
tx=fv(1);
%sl=length(spektra);
ty=min(spektrd);
text(tx,ty*3,['M_0= ',mostr,' Nm;','  {\it f}_c= ',fcstr1,' Hz;',...
    ' {\it f}_m_a_x= ', fcstr2,' Hz;',' p = ' ,mstrp,' k= ' ,strk,],'HorizontalAlignment','left',...
    'BackgroundColor',[1 1 .9],'EdgeColor','r','LineWidth',1,'FontSize',12); 
text(tx,ty,[' M_w= ',mstr1,';','{\it r} = ',fcstrr,' m;', '\Delta\sigma=',mstr2,' MPa'],...
    'HorizontalAlignment','left','BackgroundColor',[1 1 .9],'EdgeColor','r',...
        'LineWidth',1,'FontSize',12);
	
    saveas(gcf,'fmaxdispautofit.fig');saveas(gcf,'dis_spectra.jpg');hold on;close
	
	
	
	%for k = 1:3; 
     
       
        % temp=['fig',num2str(k),'.jpg']; 
        %   saveas(gca,temp); 
       %end




end
end
end