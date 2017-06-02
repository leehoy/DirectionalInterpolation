% Directional interpolation method for CT sinogram
% Reference: Directional View Interpolation for Compensation of Sparse Angular Sampling in Cone-Beam CT,¡± _IEEE Transactions on Medical Imaging_ , vol. 28, no. 7, p. 1011, Jan. 2009.
% Author: Hoyoen Lee
% Date: Dec., 07, 2016
nx=750;
% ny=200;
ny=170;
f=fopen('Sinogram\170\NoiseFree\0088.dat');
sino=fread(f,[nx ny],'float32');
fclose(f);
M=170;
B=44.5059/2;
PHI=2*pi;
a=44.5059/nx;
d_umax=B*PHI/M;
r=d_umax/a;
sigma=5;
window=10;
% projection=zeros(nx,ny*2);
% yq=1:ny*2;
% y=1:2:ny*2;
% for i=1:nx
%     projection(i,:)=interp1(y,sino(i,:),yq,'linear','extrap');
% end
projection=sino;
[J11,J12,J21,J22]=ImageOrientation2D_conv(projection,sigma,window);
interpol_window=5;
c=0.1;
treshold1=0.001;
Kmax=8;
NewProjection=zeros(nx,ny*2);
NewProjection(:,1:2:end)=projection(:,:);
RealCoord_u=-B:a:B;
RealCoord_w=linspace(0,2*pi,ny+1);
RealCoord_u=RealCoord_u(1:end-1);
RealCoord_w=RealCoord_w(1:end-1);
RealCoord_w=RealCoord_w*(1.01/(2*pi));
% for k=1:nz
%     NewProjection(:,:,(k-1)*2+1)=projection(:,:,k);
% end
for j=2:2:ny*2
    for i=1:nx
        x=(i-nx/2)*d_umax;
        y=(j/2+0.5)*B/(2*pi);
        InterpolX=i-interpol_window+1:i+interpol_window;
        InterpolY=(j/2)-interpol_window+1:(j/2)+interpol_window;
        if(j==ny*2)
            J=[(J11(i,j/2)+J11(i,1))/2,(J12(i,j/2)+J12(i,1))/2;(J21(i,j/2)+J21(i,1))/2,(J22(i,j/2)+J22(i,1))/2];
        else
            J=[(J11(i,j/2)+J11(i,j/2+1))/2,(J12(i,j/2)+J12(i,j/2+1))/2;(J21(i,j/2)+J21(i,j/2+1))/2,(J22(i,j/2)+J22(i,j/2+1))/2];
        end
        [vec,val]=eig(J);
%         InterpolX(InterpolX>nx)=InterpolX(InterpolX>nx)-nx;
%         InterpolX(InterpolX<1)=InterpolX(InterpolX<1)+nx;
%         InterpolY(InterpolY<1)=InterpolY(InterpolY<1)+ny;
%         InterpolY(InterpolY>ny)=InterpolY(InterpolY>ny)-ny;
        [xx,yy]=meshgrid((InterpolX-nx/2)*d_umax,InterpolY*B/(2*pi));
        xx=xx-x;yy=yy-y;
        d=sqrt(xx.^2+yy.^2);
        xx(d==0)=0;
        yy(d==0)=0;
        d(d==0)=1;
        xx=xx./d;yy=yy./d; % make the vectors to unit vector
        InterpolVectors=vec; % eigen vector of a pixel
        InterpolValues=val; %eigen value of a pixel
        if(max(InterpolValues(:))<treshold1)
            y_point=InterpolY;
            InterpolY(InterpolY<1)=InterpolY(InterpolY<1)+ny;
            InterpolY(InterpolY>ny)=InterpolY(InterpolY>ny)-ny;
            vq=interp1(y_point,projection(i,InterpolY),j);
            NewProjection(i,j)=vq;
            continue;
        end
        InterpolX(InterpolX>nx)=InterpolX(InterpolX>nx)-nx;
        InterpolX(InterpolX<1)=InterpolX(InterpolX<1)+nx;
        InterpolY(InterpolY<1)=InterpolY(InterpolY<1)+ny;
        InterpolY(InterpolY>ny)=InterpolY(InterpolY>ny)-ny;
        nn=zeros(size(xx));
        dd=zeros(size(xx));
        for kk=1:2
            tmp=zeros(size(xx));
            tmp=tmp+InterpolValues(kk,kk).*InterpolVectors(1,kk).*xx;
            tmp=tmp+InterpolValues(kk,kk).*InterpolVectors(2,kk).*yy;%inner product of a vector
%             tmp=tmp+InterpolValues(:,:,kk,kk).*InterpolVectors(:,:,kk,3).*zz; 
            nn=nn+tmp.^2; %sum them up for numerator
            dd=dd+InterpolValues(kk,kk).^2; % summation of eigen value for denominaor
        end
        
        dd=dd*c;
        nn(dd==0)=0;
        dd(dd==0)=1;
        weights=((1+nn./dd).^(-2))./d;
        aa=weights.*projection(InterpolX,InterpolY);
        NewProjection(i,j)=sum(aa(:))/sum(weights(:));  
    end
end
imshow(NewProjection,[]);

ny=size(NewProjection,2);
projection2=zeros(nx,ny*2);
yq=1:ny*2;
y=1:2:2*ny;
% for i=1:nx
%     projection2(i,:)=interp1(y,NewProjection(i,:),yq,'linear','extrap');
% end
% sigma=5;
% window=10;
% nz=size(projection2,3);
projection2=NewProjection;
[J11,J12,J21,J22]=ImageOrientation2D_conv(projection2,sigma,window);
interpol_window=5;
% c=0.5;
% treshold1=0.001;
% Kmax=8;
NewProjection_Dir=zeros(nx,ny*2);
NewProjection_Dir(:,1:2:end)=projection2(:,:);
for j=2:2:ny*2
    for i=1:nx
        x=i;
        y=j;
        InterpolX=i-interpol_window+1:i+interpol_window;
        InterpolY=(j)-interpol_window+1:(j)+interpol_window;
        J=[J11(i,j),J12(i,j);J21(i,j),J22(i,j)];
        [vec,val]=eig(J);
        [xx,yy]=meshgrid(InterpolX,InterpolY);
        xx=xx-x;yy=yy-y;
        d=sqrt(xx.^2+yy.^2);
        xx(d==0)=0;
        yy(d==0)=0;
        d(d==0)=1;
        xx=xx./d;yy=yy./d; % make the vectors to unit vector
        InterpolVectors=vec; % eigen vector of a pixel
        InterpolValues=val; %eigen value of a pixel
        if(max(InterpolValues(:))<treshold1)
            y_point=InterpolY;
            InterpolY(InterpolY<1)=InterpolY(InterpolY<1)+ny;
            InterpolY(InterpolY>ny)=InterpolY(InterpolY>ny)-ny;
            vq=interp1(y_point,projection2(i,InterpolY),j);
            NewProjection_Dir(i,j)=vq;
            continue;
        end
        InterpolX(InterpolX>nx)=InterpolX(InterpolX>nx)-nx;
        InterpolX(InterpolX<1)=InterpolX(InterpolX<1)+nx;
        InterpolY(InterpolY<1)=InterpolY(InterpolY<1)+ny;
        InterpolY(InterpolY>ny)=InterpolY(InterpolY>ny*2)-ny;
        nn=zeros(size(xx));
        dd=zeros(size(xx));
        for kk=1:2
            tmp=zeros(size(xx));
            tmp=tmp+InterpolValues(kk,kk).*InterpolVectors(1,kk).*xx;
            tmp=tmp+InterpolValues(kk,kk).*InterpolVectors(2,kk).*yy;%inner product of a vector
            nn=tmp.^2; %sum them up for numerator
            dd=dd+InterpolValues(kk,kk).^2; % summation of eigen value for denominaor
        end
        dd=dd*c;
        nn(dd==0)=0;
        dd(dd==0)=1;
        weights=((1+nn./dd).^(-2))./d;
        aa=weights.*projection2(InterpolX,InterpolY);
        NewProjection_Dir(i,j)=sum(aa(:))/sum(weights(:));  
    end
end
figure; imshow(NewProjection_Dir,[]);
% % write_path='D:\DeepLearning\DirectionalInterpolation\Interp170\'; %Write path
% write_path='.\'; %Write path
% sino=zeros(nx,680,ny);
% for i=1:ny
%     sino(:,:,i)=NewProjection(:,i,:);
% end
% for i=1:ny
%     f=fopen(strcat(write_path,sprintf('%04d.dat',i-1)),'w');
%     fwrite(f,sino(:,:,i),'float32');
%     fclose(f);
% end