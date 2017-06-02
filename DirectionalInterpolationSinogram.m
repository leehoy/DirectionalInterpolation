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
sino=projection;
% path='D:\DeepLearning\Sino170\'; %Sinogram directory
% list=dir(strcat(path,'*.dat'));
% sino=zeros(nx,nz,ny);
% for i=1:length(list)
%     f=fopen(strcat(path,list(i).name));
%     sino(:,:,i)=fread(f,[nx nz],'float32');
%     fclose(f);
% end
% projection=zeros(nx,ny,nz);
% for i=1:nz
%     projection(:,:,i)=sino(:,i,:);
% end
% imshow(sino(:,:,10),[])
sigma=5;
window=10;
[Vectors,Values]=ImageOrientation2D(projection,sigma,window);
interpol_window=5;
c=0.5;
treshold1=0.001;
Kmax=8;
NewProjection=zeros(nx,ny*2);
NewProjection(:,1:2:end)=projection(:,:);
% for k=1:nz
%     NewProjection(:,:,(k-1)*2+1)=projection(:,:,k);
% end
for j=2:2:ny*2
    for i=1:nx
        x=i;
        y=j/2+0.5;
        InterpolX=i-interpol_window+1:i+interpol_window;
        InterpolY=(j/2)-interpol_window+1:(j/2)+interpol_window;
%         InterpolX(InterpolX>nx)=InterpolX(InterpolX>nx)-nx;
%         InterpolX(InterpolX<1)=InterpolX(InterpolX<1)+nx;
%         InterpolY(InterpolY<1)=InterpolY(InterpolY<1)+ny;
%         InterpolY(InterpolY>ny)=InterpolY(InterpolY>ny)-ny;
        [xx,yy]=meshgrid(InterpolX,InterpolY);
        xx=xx-x;yy=yy-y;
        d=sqrt(xx.^2+yy.^2);
        xx=xx./d;yy=yy./d; % make the vectors to unit vector
        InterpolVectors=Vectors(i,j/2,:,:); % eigen vector of a pixel
        InterpolValues=Values(i,j/2,:,:); %eigen value of a pixel
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
            tmp=tmp+InterpolValues(1,1,kk,kk).*InterpolVectors(1,1,1,kk).*xx;
            tmp=tmp+InterpolValues(1,1,kk,kk).*InterpolVectors(1,1,2,kk).*yy;%inner product of a vector
%             tmp=tmp+InterpolValues(:,:,kk,kk).*InterpolVectors(:,:,kk,3).*zz; 
            nn=tmp.^2; %sum them up for numerator
            dd=dd+InterpolValues(1,1,kk,kk).^2; % summation of eigen value for denominaor
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
[Vectors,Values]=ImageOrientation(projection2,sigma,window);
interpol_window=5;
% c=0.5;
% treshold1=0.001;
% Kmax=8;
NewProjection=zeros(nx,ny*2);
NewProjection(:,1:2:end)=projection2(:,:);
for j=2:2:ny*2
    for i=1:nx
        x=i;
        y=j/2+0.5;
        InterpolX=i-interpol_window+1:i+interpol_window;
        InterpolY=(j/2)-interpol_window+1:(j/2)+interpol_window;
        [xx,yy]=meshgrid(InterpolX,InterpolY);
        xx=xx-x;yy=yy-y;
        d=sqrt(xx.^2+yy.^2);
        xx=xx./d;yy=yy./d; % make the vectors to unit vector
        InterpolVectors=Vectors(i,j/2,:,:); % eigen vector of a pixel
        InterpolValues=Values(i,j/2,:,:); %eigen value of a pixel
        if(max(InterpolValues(:))<treshold1)
            y_point=InterpolY;
            InterpolY(InterpolY<1)=InterpolY(InterpolY<1)+ny;
            InterpolY(InterpolY>ny)=InterpolY(InterpolY>ny)-ny;
            vq=interp1(y_point,projection2(x,InterpolY),y);
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
            tmp=tmp+InterpolValues(1,1,kk,kk).*InterpolVectors(1,1,1,kk).*xx;
            tmp=tmp+InterpolValues(1,1,kk,kk).*InterpolVectors(1,1,2,kk).*yy;%inner product of a vector
            nn=tmp.^2; %sum them up for numerator
            dd=dd+InterpolValues(1,1,kk,kk).^2; % summation of eigen value for denominaor
        end
        dd=dd*c;
        nn(dd==0)=0;
        dd(dd==0)=1;
        weights=((1+nn./dd).^(-2))./d;
        aa=weights.*projection2(InterpolX,InterpolY);
        NewProjection(i,j)=sum(aa(:))/sum(weights(:));  
    end
end
imshow(NewProjection,[]);
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