function []=Optic_Chiasm;

%========================
%Geometrical definition
%========================

%Voxels and Downsampling
steps=256;
bin=4;
bins=ceil(steps/bin);

SNR=20;

%Generalized sigmoid function parameters
A=0.05;       %Curve scaling factor
T=10;

Centre=[128,128];

%Set transition curvature where sigmoid fibre should transition to curved
dY=0.025;

%dY=(A*exp(-A*x)/(1+A*exp(-A*x))^2);
%U=((2-1/dY)+sqrt(1/dY^2 -4/dY))/2;
%XU=real((log(A/U))/A);

%XU=real((log(A/(((2-1/dY)+sqrt(1/dY^2 -4/dY))/2)))/A);

%fx=1/(1+exp(-n*x);
%dfx=-n/(1+exp(-n*x) * (1/(1+exp(-n*x) -1);
%Dsq=(i-x)^2 + (j-(1/(1+exp(-A*x))))^2;



%Lesion stats unused
y_les=100;
x_les=150;


%========================
% Biophysical system
%========================

%fibres details  (%Gaussian for now but can be be anything)
LAMBDA=1.7;
BETA  =0.2;

Dcsf=2.5;
Dex=0.7;

%fiber volume fractions are pre defined
fr1=0.5;
fr2=0.5;
fr3=0.5;

%========================
% MRI system
%========================

X=ndirez(60);

for i=1:60
    X(i,:)=X(i,:)/norm(X(i,:));
end

%bvals
bval(1:66)=3000;   
bval(1:6)=0;    
bvals=bval;
%bvecs
bvec(1:66,1:3)=0;
bvec(7:end,:)=X;
bvecs=bvec';

vols=length(bval);

DATA=zeros(steps,steps,10,vols);
FF=zeros(steps,steps,10,5);

DATAs=zeros(bins+2,bins+2,10,vols);
FFs=zeros(bins+2,bins+2,10,5);

f1=zeros(steps,steps);
f2=zeros(steps,steps);
f3=zeros(steps,steps);
f4=zeros(steps,steps);

S1=zeros(steps,steps,vols);
S2=zeros(steps,steps,vols);
S3=zeros(steps,steps,vols);
S4=zeros(steps,steps,vols);

%============================
%signal generation
%============================

for i=1:steps
    for j=1:steps
        
        %reverse direction
        k=(steps-j);
            
        
    	%Circle A
%        if  ((steps-j)-OA_X)^2 / (a+T)^2 + (i-OA_Y)^2 / (b+T)^2 < 1 && ((steps-j)-OA_X)^2 / a^2 + (i-OA_Y)^2 / b^2 > 1
%            
%            vect1=[(steps-j)-OA_X -(i-OA_Y)]/norm([(steps-j)-OA_X -(i-OA_Y)]);
%            fibdir1=[-vect1(2) vect1(1) 0];
%            
%            S1(j,i,1:6)=1;
%            S1(j,i,7:66)=sign_gen(fibdir1,X,3,LAMBDA,BETA);
%            f1(j,i)=fr1;
%        end
          
        %Circle B
%        if ((steps-j)-OB_X)^2 / (a2+T2)^2 + (i-OB_Y)^2 / (b2+T2)^2 < 1 && ((steps-j)-OB_X)^2 / a2^2 + (i-OB_Y)^2 / b2^2 > 1
%            
%            vect2=[(steps-j)-OB_X -(i-OB_Y)]/norm([(steps-j)-OB_X -(i-OB_Y)]);
%            fibdir2=[vect2(2) -vect2(1) 0];
%            
%            S2(j,i,1:6)=1;
%            S2(j,i,7:66)=sign_gen(fibdir2,X,3,LAMBDA,BETA);
%            f2(j,i)=fr2;
%        end
       
        %Sigmoid is centred around 0, we shift it 128 units right
        %Shift up by Thickness to avoid clipping
        k2=k-128;
        i2=i-T;
        
        
        %Functions for sigmoid curves
        y=@(x)(steps-2*T)/(1+exp(-A*x));
        y2=@(x)(steps-2*T)/(1+exp(A*x));
        
        %Dsq is distance function of a point i, k2 from the curve
        Dsq=@(x)((k2-x)^2) + (i2-y(x))^2;
        Dsq2=@(x)((k2-x)^2) + (i2-y2(x))^2;
        
        %p0 finds closest x on curve from point being considered
        p0=fminsearch(Dsq,k2);
        p02=fminsearch(Dsq2,k2);
        
        if (k2-p0)^2 + (i2-y(p0))^2 < T^2
            grad = (steps-2*T)*(A*exp(-A*p0))/((1+exp(-A*p0))^2);
            vect1=[1 -grad];
            vect1=vect1/norm(vect1);
            
            fibdir1=[vect1(1) vect1(2) 0];
            
            S1(j,i,1:6)=1;
            S1(j,i,7:66)=sign_gen(fibdir1,X,3,LAMBDA,BETA);
            f1(j,i)=fr1;
        end
        
        if (k2-p02)^2 + (i2-y2(p02))^2 < T^2
            grad=(steps-2*T)*(A*exp(A*p02))/((1+exp(A*p02))^2);
            vect2=[1 -grad];
            vect2=vect2/norm(vect2);
            fibdir2=[vect2(1) -vect2(2) 0];
            
            S2(j,i,1:6)=1;
            S2(j,i,7:66)=sign_gen(fibdir2,X,3,LAMBDA,BETA);
            f2(j,i)=fr2;
        end       
        
        if f1(j,i)+f2(j,i)+f3(j,i)>0 && f1(j,i)+f2(j,i)+f3(j,i)<1
            f4(j,i)= 1 - (f1(j,i)+f2(j,i)+f3(j,i));       
            S4(j,i,1:6)=1;
            S4(j,i,7:66)=sign_gen([1 1 0],X,3,Dcsf,Dcsf);
        end
    end
end

%%%%%%%%%%
%Combining
%%%%%%%%%%
for t=3:8
    FF(:,:,t,1)=f1;
    FF(:,:,t,2)=f2;
    FF(:,:,t,3)=f3;  
    FF(:,:,t,4)=f4;
    
    DATA(:,:,t,:)=(f1.*S1+f2.*S2+f3.*S3+f4.*S4);
end

for i=1:length(DATA(1,1,:,1))
    for j=1:length(DATA(1,1,1,:))
      DATAs(:,:,i,j)=sgn_noise(DATAs(:,:,i,j),SNR);
    end
end

%Gradient at centre where p0=0 so exp(-Ax) = 1
%Angle between Fibres given by 180 - 2atan(grad)
CentreGrad = (steps-2*T)*A/4;
Angle = pi-2*atan(CentreGrad);
AngleDeg=Angle*180/pi;

data_name = sprintf('Optic_Chiasm%.2f_%dSNR_Data.nii',AngleDeg,SNR);
VF_name = sprintf('Optic_Chiasm%.2f_%dSNR_VolumeFractions.nii',AngleDeg,SNR);
mask_name = sprintf('Optic_Chiasm%.2f_%dSNR_DataMask.nii',AngleDeg,SNR);
bvals_name = sprintf('Optic_Chiasm%.2f_%dSNR.bvals',AngleDeg,SNR);
bvecs_name = sprintf('Optic_Chiasm%.2f_%dSNR.bvecs',AngleDeg,SNR);

%%%%%%%%%%
%size down
%%%%%%%%%%

DATAs(2:end-1,2:end-1,:,:)=imresize(DATA,1/bin,'bilinear');
FFs(2:end-1,2:end-1,:,:)=imresize(FF,1/bin,'bilinear');


%=================================
%saving data at high resolution
%=================================

DATAnii=make_nii(single(DATAs),[2 2 2]);
save_nii(DATAnii,data_name);

FFnii=make_nii(single(FFs),[2 2 2]);
save_nii(FFnii,VF_name);

mask=DATAs(:,:,:,1);
masknii=make_nii(single(mask),[2 2 2]);
save_nii(masknii,mask_name);

save(bvals_name,'bvals','-ASCII');
save(bvecs_name,'bvecs','-ASCII');








