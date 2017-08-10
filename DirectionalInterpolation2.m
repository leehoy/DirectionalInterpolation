% Directional interpolation method for CT sinogram
% Reference: Directional View Interpolation for Compensation of Sparse Angular Sampling in Cone-Beam CT,¡± _IEEE Transactions on Medical Imaging_ , vol. 28, no. 7, p. 1011, Jan. 2009.
% Author: Hoyoen Lee
% Date: Dec., 07, 2016
nx=750;
ny=200;
nz=170;
path='D:\DeepLearning\Sino170\'; %Sinogram directory
list=dir(strcat(path,'*.dat'));
sino=zeros(nx,nz,ny);
for i=1:length(list)
    f=fopen(strcat(path,list(i).name));
    sino(:,:,i)=fread(f,[nx nz],'float32');
    fclose(f);
end
projection=zeros(nx,ny,nz);
for i=1:nz
    projection(:,:,i)=sino(:,i,:);
end
imshow(sino(:,:,10),[])
sigma=1;
window=10;
[Vectors,Values]=ImageOrientation(projection,sigma,window);
interpol_window=5;
c=0.5;
treshold1=0.001;

NewProjection=zeros(nx,ny,nz*2);
for k=1:nz
    NewProjection(:,:,(k-1)*2+1)=projection(:,:,k);
end
for k=2:2:nz*2
    for i=1:nx
        for j=1:ny
            InterpolX=i-interpol_window+1:i+interpol_window;
            InterpolY=j-interpol_window+1:j+interpol_window;
            InterpolZ=[k-1,k+1];
            InterpolX(InterpolX>nx)=nx; %InterpolX(InterpolX>nx)-nx;
            InterpolX(InterpolX<1)=1; %InterpolX(InterpolX<1)+nx;
            InterpolY(InterpolY<1)=1; %InterpolY(InterpolY<1)+ny;
            InterpolY(InterpolY>ny)=ny; %InterpolY(InterpolY>ny)-ny;
            InterpolZ(InterpolZ<1)=InterpolZ(InterpolZ<1)+nz;
            InterpolZ(InterpolZ>(nz*2))=InterpolZ(InterpolZ>(nz*2))-(nz*2);
            
            [xx,yy,zz]=meshgrid(InterpolX,InterpolY,InterpolZ);
            xx=xx-i;yy=yy-j;zz=zz-k;
            d=sqrt(xx.^2+yy.^2+zz.^2);
            xx=xx./d;yy=yy./d;zz=zz./d;
            InterpolVectors=Vectors(InterpolX,InterpolY,InterpolZ,:,:);
            InterpolValues=Values(InterpolX,InterpolY,InterpolZ,:,:);
            if(max(InterpolValues(:))<treshold1)
                continue;
            end
            nn=zeros(size(xx));
            dd=zeros(size(xx));
            for kk=1:3
                tmp=zeros(size(xx));
                tmp=tmp+InterpolValues(:,:,:,kk,kk).*InterpolVectors(:,:,:,kk,1).*xx;
                tmp=tmp+InterpolValues(:,:,:,kk,kk).*InterpolVectors(:,:,:,kk,2).*yy;
                tmp=tmp+InterpolValues(:,:,:,kk,kk).*InterpolVectors(:,:,:,kk,3).*zz;
                nn=nn+tmp.^2;
                dd=dd+InterpolValues(:,:,:,kk,kk).^2;
            end
            dd=dd*c;
            nn(dd==0)=0;
            dd(dd==0)=1;
            weights=((1+nn./dd).^(-2))./d;
            aa=weights.*NewProjection(InterpolX,InterpolY,InterpolZ);
            NewProjection(i,j,k)=sum(aa(:))/sum(weights(:));  
        end
    end
end
imshow(NewProjection(:,:,4),[]);
projection2=NewProjection;
sigma=1;
window=10;
nz=size(projection2,3);
[Vectors,Values]=ImageOrientation(projection2,sigma,window);
interpol_window=5;
c=0.5;
treshold1=0.001;
Kmax=8;
NewProjection=zeros(nx,ny,nz*2);
for k=1:nz
    NewProjection(:,:,(k-1)*2+1)=projection2(:,:,k);
end
for k=2:2:nz*2
    for i=1:nx
        for j=1:ny
            InterpolX=i-interpol_window+1:i+interpol_window;
            InterpolY=j-interpol_window+1:j+interpol_window;
            InterpolZ=[k-1,k+1];
            InterpolX(InterpolX>nx)=nx; %InterpolX(InterpolX>nx)-nx;
            InterpolX(InterpolX<1)=1; %InterpolX(InterpolX<1)+nx;
            InterpolY(InterpolY<1)=1; %InterpolY(InterpolY<1)+ny;
            InterpolY(InterpolY>ny)=ny; %InterpolY(InterpolY>ny)-ny;
            InterpolZ(InterpolZ<1)=InterpolZ(InterpolZ<1)+nz;
            InterpolZ(InterpolZ>(nz*2))=InterpolZ(InterpolZ>(nz*2))-(nz*2);
            
            [xx,yy,zz]=meshgrid(InterpolX,InterpolY,InterpolZ);
            xx=xx-i;yy=yy-j;zz=zz-k;
            d=sqrt(xx.^2+yy.^2+zz.^2);
            xx=xx./d;yy=yy./d;zz=zz./d;
            InterpolVectors=Vectors(InterpolX,InterpolY,InterpolZ,:,:);
            InterpolValues=Values(InterpolX,InterpolY,InterpolZ,:,:);
            if(max(InterpolValues(:))<treshold1)
                continue;
            end
            nn=zeros(size(xx));
            dd=zeros(size(xx));
            for kk=1:3
                tmp=zeros(size(xx));
                tmp=tmp+InterpolValues(:,:,:,kk,kk).*InterpolVectors(:,:,:,kk,1).*xx;
                tmp=tmp+InterpolValues(:,:,:,kk,kk).*InterpolVectors(:,:,:,kk,2).*yy;
                tmp=tmp+InterpolValues(:,:,:,kk,kk).*InterpolVectors(:,:,:,kk,3).*zz;
                nn=nn+tmp.^2;
                dd=dd+InterpolValues(:,:,:,kk,kk).^2;
            end
            dd=dd*c;
            nn(dd==0)=0;
            dd(dd==0)=1;
            weights=((1+nn./dd).^(-2))./d;
            aa=weights.*NewProjection(InterpolX,InterpolY,InterpolZ);
            NewProjection(i,j,k)=sum(aa(:))/sum(weights(:));  
        end
    end
end
imshow(NewProjection(:,:,4),[]);
% write_path='D:\DeepLearning\DirectionalInterpolation\Interp170\'; %Write path
write_path='.\'; %Write path
sino=zeros(nx,680,ny);
for i=1:ny
    sino(:,:,i)=NewProjection(:,i,:);
end
for i=1:ny
    f=fopen(strcat(write_path,sprintf('%04d.dat',i-1)),'w');
    fwrite(f,sino(:,:,i),'float32');
    fclose(f);
end