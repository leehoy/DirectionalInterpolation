% Directional interpolation method for CT sinogram
% Reference: Directional View Interpolation for Compensation of Sparse Angular Sampling in Cone-Beam CT,¡± _IEEE Transactions on Medical Imaging_ , vol. 28, no. 7, p. 1011, Jan. 2009.
% Author: Hoyoen Lee
% Date: Dec., 07, 2016
tic;
nx=750;
% ny=200;
ny=170;
f=fopen('Sinogram\170\NoiseFree\0088.dat');
projection=fread(f,[nx ny],'float32');
fclose(f);

sigma=0.05;
window=20;
% [Vectors,Values]=ImageOrientation2D(projection,sigma,window);
[J11,J12,J21,J22]=ImageOrientation2D_conv(projection,sigma,window);
interpol_window=5;
c=0.01;
treshold1=0.001;
Kmax=8;
NewProjection=zeros(nx,ny*2);
NewProjection(:,1:2:end)=projection(:,:);
ys=45.0/(2*pi);
xs=44.5/750;
for i=1:nx
    for j=1:ny
        x=i;
        y=2*j;
        InterpolX=i-interpol_window:i+interpol_window;
        InterpolY=j-interpol_window:j+interpol_window;
        [xx,yy]=meshgrid(InterpolX,InterpolY*ys);
        xx=xx*xs-x;yy=yy-(j+0.5)*ys;
        d=sqrt(xx.^2+yy.^2);
        d(d==0)=1;
        xx=xx./d;yy=yy./d; % make the vectors to unit vector
        if(j>=ny)
            J=[(J11(i,1)+J11(i,j))/2,(J12(i,1)+J12(i,j))/2;(J21(i,1)+J21(i,j))/2,(J22(i,1)+J22(i,j))/2];
        else
            J=[(J11(i,j)+J11(i,j+1))/2,(J12(i,j)+J12(i,j+1))/2;(J21(i,j)+J21(i,j+1))/2,(J22(i,j)+J22(i,j+1))/2];
        end
        [V,D]=eig(J);
        InterpolVectors=V;
        InterpolValues=D;
        if(max(InterpolValues(:))<treshold1)
            y_point=InterpolY*ys;
            InterpolY(InterpolY<1)=InterpolY(InterpolY<1)+ny;
            InterpolY(InterpolY>ny)=InterpolY(InterpolY>ny)-ny;
            vq=interp1(y_point,projection(i,InterpolY),(j+0.5)*ys,'linear','extrap');
            NewProjection(x,y)=vq;
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
            nn=nn+tmp.^2; %sum them up for numerator
            dd=dd+InterpolValues(kk,kk).^2; % summation of eigen value for denominaor
        end
        
        dd=dd*c;
        nn(dd==0)=0;
        dd(dd==0)=1;
        weights=((1+nn./dd).^(-2))./d;
        aa=weights.*projection(InterpolX,InterpolY);
        NewProjection(x,y)=sum(aa(:))/sum(weights(:));  
    end
end
close all;
imshow(NewProjection,[]);
toc

nx=750;
ny=340;
projection=NewProjection;

sigma=0.05;
window=20;
% [Vectors,Values]=ImageOrientation2D(projection,sigma,window);
[J11,J12,J21,J22]=ImageOrientation2D_conv(projection,sigma,window);
interpol_window=5;
c=0.01;
treshold1=0.001;
Kmax=8;
NewProjection2=zeros(nx,ny*2);
NewProjection2(:,1:2:end)=projection(:,:);
ys=45.0/(2*pi);
xs=44.5/750;
for i=1:nx
    for j=1:ny
        x=i;
        y=2*j;
        InterpolX=i-interpol_window:i+interpol_window;
        InterpolY=j-interpol_window:j+interpol_window;
        [xx,yy]=meshgrid(InterpolX,InterpolY*ys);
        xx=xx*xs-x;yy=yy-(j+0.5)*ys;
        d=sqrt(xx.^2+yy.^2);
        d(d==0)=1;
        xx=xx./d;yy=yy./d; % make the vectors to unit vector
        if(j>=ny)
            J=[(J11(i,1)+J11(i,j))/2,(J12(i,1)+J12(i,j))/2;(J21(i,1)+J21(i,j))/2,(J22(i,1)+J22(i,j))/2];
        else
            J=[(J11(i,j)+J11(i,j+1))/2,(J12(i,j)+J12(i,j+1))/2;(J21(i,j)+J21(i,j+1))/2,(J22(i,j)+J22(i,j+1))/2];
        end
        [V,D]=eig(J);
        InterpolVectors=V;
        InterpolValues=D;
        if(max(InterpolValues(:))<treshold1)
            y_point=InterpolY*ys;
            InterpolY(InterpolY<1)=InterpolY(InterpolY<1)+ny;
            InterpolY(InterpolY>ny)=InterpolY(InterpolY>ny)-ny;
            vq=interp1(y_point,projection(i,InterpolY),(j+0.5)*ys,'linear','extrap');
            NewProjection2(x,y)=vq;
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
            nn=nn+tmp.^2; %sum them up for numerator
            dd=dd+InterpolValues(kk,kk).^2; % summation of eigen value for denominaor
        end
        
        dd=dd*c;
        nn(dd==0)=0;
        dd(dd==0)=1;
        weights=((1+nn./dd).^(-2))./d;
        aa=weights.*projection(InterpolX,InterpolY);
        NewProjection2(x,y)=sum(aa(:))/sum(weights(:));  
    end
end
close all;
imshow(NewProjection2,[])