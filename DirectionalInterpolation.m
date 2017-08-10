function [ NewProjection ] = DirectionalInterpolation( image,NumberOfUpscale,sigma_x,sigma_y,interpol_window,scale)
%UNTITLED2 이 함수의 요약 설명 위치
%   자세한 설명 위치

c=0.01;
treshold=0.001;
xs=1;
for k=1:NumberOfUpscale
    nx=size(image,1);
    ny=size(image,2);
    NewProjection=zeros(nx,ny*2);
    NewProjection(:,1:2:end)=image(:,:);
    ys=scale/(2*pi);
    [J11,J12,J21,J22]=ImageOrientation2D_conv(image,sigma_x,sigma_y);
    for i=1:nx
        for j=1:ny
            x=i;
            y=2*j;
            InterpolX=i-interpol_window:i+interpol_window;
            InterpolY=j:j+1;
            [yy,xx]=meshgrid(InterpolY*ys,InterpolX*xs);
            xx=xx-x*xs;yy=yy-(j+0.5)*ys;
            d=sqrt(xx.^2+yy.^2);
            d(d==0)=1;
            xx=xx./d;yy=yy./d; % make the vectors to unit vector
            if(j>=ny)
                J=[(J11(i,1)+J11(i,j))/2,(J12(i,1)+J12(i,j))/2;(J21(i,1)+J21(i,j))/2,(J22(i,1)+J22(i,j))/2];
            else
                J=[(J11(i,j)*0.5*ys+J11(i,j+1)*0.5*ys)/ys,(J12(i,j)*0.5*ys+J12(i,j+1)*0.5*ys)/ys;(J21(i,j)*0.5*ys+J21(i,j+1)*0.5*ys)/ys,(J22(i,j)*0.5*ys+J22(i,j+1)*0.5*ys)/ys];
            end
            [V,D]=eig(J);
            InterpolVectors=V;
            InterpolValues=D;
            if(max(InterpolValues(:))<treshold)
                y_point=InterpolY*ys;
                InterpolY(InterpolY<1)=InterpolY(InterpolY<1)+ny;
                InterpolY(InterpolY>ny)=InterpolY(InterpolY>ny)-ny;
                vq=interp1(y_point,image(i,InterpolY),(j+0.5)*ys,'linear','extrap');
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
            aa=weights.*image(InterpolX,InterpolY);
            NewProjection(x,y)=sum(aa(:))/sum(weights(:));  
        end
    end
    image=NewProjection;
end
end

