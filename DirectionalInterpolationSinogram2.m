% Directional interpolation method for CT sinogram
% Reference: Directional View Interpolation for Compensation of Sparse Angular Sampling in Cone-Beam CT,¡± _IEEE Transactions on Medical Imaging_ , vol. 28, no. 7, p. 1011, Jan. 2009.
% Author: Hoyoen Lee
% Date: Dec., 07, 2016
nx=750;
% ny=200;
ny=170;
f=fopen('Sinogram\170\NoiseFree\0088.dat');
projection=fread(f,[nx ny],'float32');
fclose(f);
sino=projection';
M=170;
B=44.5059/2;
PHI=2*pi;
a=44.5059/nx;
d_umax=B*PHI/M;
r=d_umax/a;
sigma=5;
window=10;
[J11,J12,J21,J22]=ImageOrientation2D_conv(projection,sigma,window);
interpol_window=5;
c=0.5;
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
        x=i;
        y=j/2+0.5;
        InterpolX=i-interpol_window+1:i+interpol_window;
        InterpolY=(j/2)-interpol_window+1:(j/2)+interpol_window;
        J=[J11(i,j/2),J12(i,j/2);J21(i,j/2),J22(i,j/2)];
        [vec,val]=eig(J);
%         InterpolX(InterpolX>nx)=InterpolX(InterpolX>nx)-nx;
%         InterpolX(InterpolX<1)=InterpolX(InterpolX<1)+nx;
%         InterpolY(InterpolY<1)=InterpolY(InterpolY<1)+ny;
%         InterpolY(InterpolY>ny)=InterpolY(InterpolY>ny)-ny;
        [xx,yy]=meshgrid(InterpolX,InterpolY);
        xx=xx-x;yy=yy-y;
        d=sqrt(xx.^2+yy.^2);
        xx=xx./d;yy=yy./d; % make the vectors to unit vector
        InterpolVectors=vec; % eigen vector of a pixel
        InterpolValues=val; %eigen value of a pixel
        if(max(InterpolValues(:))<treshold1)
            y_point=InterpolY;
            InterpolY(InterpolY<1)=InterpolY(InterpolY<1)+ny;
            InterpolY(InterpolY>ny)=InterpolY(InterpolY>ny)-ny;
            vq=interp1(y_point,projection(x,InterpolY),y);
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
            nn=tmp.^2; %sum them up for numerator
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
projection2=NewProjection;
ny=size(projection2,2);
% sigma=5;
% window=10;
% nz=size(projection2,3);
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
        y=j/2+0.5;
        InterpolX=i-interpol_window+1:i+interpol_window;
        InterpolY=(j/2)-interpol_window+1:(j/2)+interpol_window;
        J=[J11(i,j/2),J12(i,j/2);J21(i,j/2),J22(i,j/2)];
        [vec,val]=eig(J);
        [xx,yy]=meshgrid(InterpolX,InterpolY);
        xx=xx-x;yy=yy-y;
        d=sqrt(xx.^2+yy.^2);
        xx=xx./d;yy=yy./d; % make the vectors to unit vector
        InterpolVectors=vec; % eigen vector of a pixel
        InterpolValues=val; %eigen value of a pixel
        if(max(InterpolValues(:))<treshold1)
            y_point=InterpolY;
            InterpolY(InterpolY<1)=InterpolY(InterpolY<1)+ny;
            InterpolY(InterpolY>ny)=InterpolY(InterpolY>ny)-ny;
            vq=interp1(y_point,projection2(x,InterpolY),y);
            NewProjection_Dir(i,j)=vq;
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