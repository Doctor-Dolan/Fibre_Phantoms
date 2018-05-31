function []=ushaped3;

%========================
%Geometrical definition
%========================


%Voxels and Downsampling
steps=256;
bin=8;
bins=ceil(steps/bin);


%=====
%Setup
%=====
%U-shaped 
O=[128 0];

a1 = 60;    %x semi axis = 1/2(a1+O(1))
b1 = 50;   

a2 = 90;
b2 = 80;

base=30;    %width of extended fibres
peak=240;   %what Y value should extended fibres reach at x=256?

Q1=O(1)-a2-base; %Y=P(x-Q)^0.5%
P1= peak/(256-Q1)^0.5;            %P calculated according to peak for Y=P(x-Q)^0.5

Q2=O(1)+a2+base;
P2= P1;                            %fibre gemoetry symmetrical for now

Q3=O(1)-a2-base/2;

SNR=0;


%theta=30;  
%thetaR=theta*pi/180;
%tan(thetaR);            %=grad
S=25;        %distance from centre for fibres to go straight

SX1=O(1)-S;
SX2=O(1)+S;

SY1=sqrt(b2^2*(1-(SX1-O(1))^2/a2^2))+O(2);
SY2=sqrt(b2^2*(1-(SX2-O(1))^2/a2^2))+O(2);        %Y coordinates where fibre goes straight

df1 = -1/((a2^2*(SY1-O(2)))/(b2^2*(SX1-O(1))));            %gradient of straight fibres 
df2 = -1/((a2^2*(SY2-O(2)))/(b2^2*(SX2-O(1))));

C1= SY1-df1*SX1;        %Y=dfX+C
C2= SY2-df2*SX2;   

%%%%%%%%%%
%filenames
%%%%%%%%%%

data_name = sprintf('UshapedExt_%d_Data.nii', SNR);
VF_name = sprintf('UshapedExt_%d_VolumeFractions.nii',SNR);
mask_name = sprintf('UshapedExt_%d_DataMask.nii',SNR);
bvals_name = sprintf('UshapedExt_%d.bvals',SNR);
bvecs_name = sprintf('UshapedExt_%d.bvecs',SNR);

%========================
% Biophysical system
%========================

%fibres details  (%Gaussian for now but can be be anything)
LAMBDA=1.7;
BETA  =0.2;

Dcsf=2.5;
Dex=0.7;



%fiber volume fractions are pre defined
fr1=0.8;
fr2=0.8;
fr3=0.8;



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
FFX=zeros(steps,steps,10,5);

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
        
        if ((steps-j)-O(1))^2/a1^2 + (i-O(2))^2/b1^2>1 && ((steps-j)-O(1))^2/a2^2 + (i-O(2))^2/b2^2<1       %U shaped ellipse
            
            ugrad=(-b2^2*((steps-j)-O(1)))/(a2^2*i);
            fibdir1=[1 -ugrad 0];
                        
            fibdir1 = fibdir1/norm(fibdir1);
            S1(j,i,1:6)=1;
            S1(j,i,7:66)=sign_gen(fibdir1,X,3,LAMBDA,BETA);
            f1(j,i)=fr1;
        end

        if ((steps-j)-O(1))^2/a2^2 + (i-O(2))^2/b2^2>1 && i<P1*((steps-j)-Q1)^0.5 && (steps-j)<=SX1     %Left hand ext
            
            bcub=4*(P1^2 -2*(steps-j)-Q1);                  %The distance function from the point ((steps-j),i) is a cubic
            ccub=(P1^4 -4*P1^2*(steps-j) + 4*(steps-j)^2 -4*P1^2*Q1 + 8*Q1*(steps-j));  %these are the coefficients (a=4)
            dcub=-Q1*(4*(steps-j)^2 + P1^4 -4*(steps-j)*P1^2) -P1^2*i^2;
            
            cubvect=double([4 bcub ccub dcub]);
            myroots=roots(cubvect);
            
            count=0;
            x_upper=0;
            for k=1:3
                if isreal(myroots(k))==1
                    count=count+1;
                    x_upper=myroots(k);
                end
            end
            if count~=1
                'houston we have a problem'     %Checking there is only one real root
            end
            
            yint=P1*(x_upper-Q1)^0.5;
            curr_dist=sqrt((yint-i)^2+(x_upper-(steps-j))^2);
            
            m=(yint-i)/(x_upper-(steps-j));       
            cc=i-m*(steps-j);                     %Line Y=mX+cc connects current point to closest point on upper boundary
            
            aquad=(m^2 + b2^2/a2^2);
            bquad=(2*m*cc -2*O(1)*b2^2/a2^2);
            cquad=(cc^2-b2^2+O(1)^2*(b2^2/a2^2));   %Solve where Y=mX+cc crosses elliptical boundary, also a polynomial.
            
            quadvect=[aquad bquad cquad];
            eroots=roots(quadvect);
            
            x_lower=0;
            count=0;
            if abs((steps-j)-eroots(1))<=abs((steps-j)-eroots(2))
                x_lower=eroots(1);
            elseif abs((steps-j)-eroots(2))<abs((steps-j)-eroots(1))
                x_lower=eroots(2);
            else
                'houston we also have a problem'
            end
            
            yeint=sqrt((1-(x_lower-O(1))^2/a2^2)*b2^2);
            
            curr_width=sqrt((yint-yeint)^2 + (x_upper-x_lower)^2);
            Qfrac=curr_dist/curr_width;                                %Current distance from curve taken normal to parabola
            
            egrad=(-b2^2*(x_lower-O(1)))/(a2^2*yeint);                  %Gradient of elliptical boundary at current point
            
            Qprime=Q1+(Qfrac*base);
            Pprime=i/((steps-j)-Qprime)^0.5;
            pgrad=Pprime/(2*((steps-j)-Qprime)^0.5);           %Gradient of parabola Y'=P'(x-Q')^0.5 Where Q' & P' are computed
                                                               %such that Qfrac is constant
                                                                        
            grad=pgrad*(1-Qfrac) + egrad*Qfrac;           %actual gradient is a weighted combination based on dist from parabola
            fibdir=[1 -grad 0];
            fibdir = fibdir/norm(fibdir);
            
            S2(j,i,1:6)=1;
            S2(j,i,7:66)=sign_gen(fibdir,X,3,LAMBDA,BETA);
            f2(j,i)=fr2*(base/curr_width);
        end    
        
        if ((steps-j)-O(1))^2/a2^2 + (i-O(2))^2/b2^2>1 && i<P2*(Q2-(steps-j))^0.5 && (steps-j)>=SX2     %Right hand ext
            
            bcub=4*(P2^2 +2*(steps-j)+Q2);                              %The distance function from the point ((steps-j),i) is a cubic
            ccub=-4*(Q2*(P2^2+2*(steps-j))+(steps-j)*(P2^2+(steps-j)))-P2^4;  %these are the coefficients (a=-4)
            dcub=Q2*(P2^4+4*(steps-j)^2+4*P2^2*(steps-j))-(i^2*P2^2);
            
            cubvect=double([-4 bcub ccub dcub]);
            myroots=roots(cubvect);
            
            count=0;
            x_upper=0;
            for k=1:3
                if isreal(myroots(k))==1
                    count=count+1;
                    x_upper=myroots(k);
                end
            end
            if count~=1
                'houston we have a problem 2'
            end
            
            yint=P2*(Q2-x_upper)^0.5;
            curr_dist=sqrt((yint-i)^2+(x_upper-(steps-j))^2);
            
            m=(yint-i)/(x_upper-(steps-j));       %Working
            cc=i-m*(steps-j);                   
            
            aquad=(m^2 + b2^2/a2^2);
            bquad=(2*m*cc -2*O(1)*b2^2/a2^2);
            cquad=(cc^2-b2^2+O(1)^2*(b2^2/a2^2));
            
            quadvect=[aquad bquad cquad];
            eroots=roots(quadvect);
            
            x_lower=0;
            count=0;
            if abs((steps-j)-eroots(1))<=abs((steps-j)-eroots(2))
                x_lower=eroots(1);
            elseif abs((steps-j)-eroots(2))<abs((steps-j)-eroots(1))
                x_lower=eroots(2);
            else
                'houston we also have a problem 3'
            end
            
            yeint=sqrt((1-(x_lower-O(1))^2/a2^2)*b2^2);
            
            curr_width=sqrt((yint-yeint)^2 + (x_upper-x_lower)^2);
            Qfrac=curr_dist/curr_width;                                %current distance from curve taken normal to parabola
            
            egrad=(-b2^2*(x_lower-O(1)))/(a2^2*yeint);
            
            Qprime=Q2-(Qfrac*base);
            Pprime=i/(Qprime-(steps-j))^0.5;
            pgrad=-Pprime/(2*(Qprime-(steps-j))^0.5);
            
            grad=pgrad*(1-Qfrac) + egrad*Qfrac;
            fibdir=[1 -grad 0];
            fibdir = fibdir/norm(fibdir);
            
            S3(j,i,1:6)=1;
            S3(j,i,7:66)=sign_gen(fibdir,X,3,LAMBDA,BETA);
            f3(j,i)=fr3*(base/curr_width);
        end 
        
        if (steps-j)>SX1 && i<P1*((steps-j)-Q1)^0.5 && i>df1*(steps-j)+C1       %Left hand ext

            bcub=4*(P1^2 -2*(steps-j)-Q1);                              %The distance function from the point ((steps-j),i) is a cubic
            ccub=(P1^4 -4*P1^2*(steps-j) + 4*(steps-j)^2 -4*P1^2*Q1 + 8*Q1*(steps-j));  %these are the coefficients (a=4)
            dcub=-Q1*(4*(steps-j)^2 + P1^4 -4*(steps-j)*P1^2) -P1^2*i^2;
            
            cubvect=double([4 bcub ccub dcub]);
            myroots=roots(cubvect);
            
            count=0;
            x_upper=0;
            for k=1:3
                if isreal(myroots(k))==1
                    count=count+1;
                    x_upper=myroots(k);
                end
            end
            if count~=1
                'houston we have a problem 7'
            end
            
            y_upper=P1*(x_upper-Q1)^0.5;
            upper_intercept=[x_upper y_upper];
            
            m=(y_upper-i)/(x_upper-(steps-j));
            cc=i-m*(steps-j);
            
            x_lower=(C1-cc)/(m-df1);
            y_lower=df1*x_lower+C1;
            lower_intercept=[x_lower y_lower];
            
            
            curr_dist=norm([(steps-j) i] - upper_intercept);
            curr_width=norm(upper_intercept-lower_intercept);
            Qfrac=curr_dist/curr_width;
            
            Qprime=Q1+(Qfrac*base);
            Pprime=i/((steps-j)-Qprime)^0.5;
            pgrad=Pprime/(2*((steps-j)-Qprime)^0.5);
            
            grad=pgrad*(1-Qfrac) + (df1*Qfrac);
            
            fibdir1 = [1 -grad 0];
            S2(j,i,1:6)=1;
            S2(j,i,7:66)=sign_gen(fibdir1,X,3,LAMBDA,BETA);
            f2(j,i)=fr2*(base/curr_width);
        end

        if (steps-j)<SX2 && i<P2*(Q2-(steps-j))^0.5 && i>df2*(steps-j)+C2       %Right hand ext
            
            bcub=4*(P2^2 +2*(steps-j)+Q2);                              %The distance function from the point ((steps-j),i) is a cubic
            ccub=-4*(Q2*(P2^2+2*(steps-j))+(steps-j)*(P2^2+(steps-j)))-P2^4;  %these are the coefficients (a=-4)
            dcub=Q2*(P2^4+4*(steps-j)^2+4*P2^2*(steps-j))-(i^2*P2^2);
            
            cubvect=double([-4 bcub ccub dcub]);
            myroots=roots(cubvect);
            
            count=0;
            x_upper=0;
            for k=1:3
                if isreal(myroots(k))==1
                    count=count+1;
                    x_upper=myroots(k);
                end
            end
            if count~=1
                'houston we have a problem 5'
            end
            
            y_upper=P2*(Q2-x_upper)^0.5;
            upper_intercept=[x_upper y_upper];
            
            m=(y_upper-i)/(x_upper-(steps-j));
            cc=i-m*(steps-j);
            
            x_lower=(C2-cc)/(m-df2);
            y_lower=df2*x_lower+C2;
            lower_intercept=[x_lower y_lower];
            
            curr_dist=norm([(steps-j) i] - upper_intercept);
            curr_width=norm(upper_intercept-lower_intercept);
            Qfrac=curr_dist/curr_width;
            
            Qprime=Q2-(Qfrac*base);
            Pprime=i/(Qprime-(steps-j))^0.5;
            pgrad=-Pprime/(2*(Qprime-(steps-j))^0.5);
            
            grad=pgrad*(1-Qfrac) + (df2*Qfrac);
            
            
            fibdir1 = [1 -grad 0];
            S3(j,i,1:6)=1;
            S3(j,i,7:66)=sign_gen(fibdir1,X,3,LAMBDA,BETA);
            f3(j,i)=fr3*(base/curr_width);
        end
        
        if f2(j,i)+f3(j,i)>0 && f2(j,i)+f3(j,i)<1
            
            f4(j,i)=1-(f1(j,i)+f2(j,i)+f3(j,i));
            S4(j,i,1:6)=1;
            S4(j,i,7:66)=sign_gen([1 1 0],X,3,Dcsf,Dcsf);
        end 
    end
end


for t=3:8
    FFX(:,:,t,1)=f1;
    FFX(:,:,t,2)=f2;
    FFX(:,:,t,3)=f3;
    FFX(:,:,t,4)=f4;
    
    DATA(:,:,t,:)=(f2.*S2);%+f3.*S3+f4.*S4);
end

%DATA=sgn_noise(DATA,SNR);

%%%%%%%%%%
%size down
%%%%%%%%%%

DATAs(2:end-1,2:end-1,:,:)=imresize(DATA,1/bin,'bilinear');
FFs(2:end-1,2:end-1,:,:)=imresize(FFX,1/bin,'bilinear');


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
