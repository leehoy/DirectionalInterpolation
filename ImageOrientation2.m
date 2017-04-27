function [ Vectors,Values] = ImageOrientation( image,sigma,window)
%Image orientation function
%   자세한 설명 위치
J=zeros([3,3]);

% [G1,G2]=imgradientxy(image(:,:,k));
% G3=image(:,:,k)-image(:,:,k+1);
Vectors=zeros(size(image,1),size(image,2),size(image,3)*2,3,3);
Values=zeros(size(image,1),size(image,2),size(image,3)*2,3,3);
PadImage=zeros(size(image,1)+2*window,size(image,2)+2*window,size(image,3));
for k=1:size(image,3)
    PadImage(:,:,k)=padarray(image(:,:,k),[window window],'circular','both');
end
x=-window+1:size(image,1)+window;
y=-window+1:size(image,2)+window;
z=-window+1:size(image,3)+window;
for k=1:size(image,3)
    [G1,G2]=imgradientxy(PadImage(:,:,k));
    if k==size(image,3)
        G3=PadImage(:,:,k)-PadImage(:,:,1);
    else
        G3=PadImage(:,:,k)-PadImage(:,:,k+1);
    end
    G11=G1.*G1;
    G12=G1.*G2;
    G13=G1.*G3;
    G21=G2.*G1;
    G22=G2.*G2;
    G23=G2.*G3;
    G31=G3.*G1;
    G32=G3.*G2;
    G33=G3.*G3;
    for i=1:size(image,1)
        for j=1:size(image,2)
            d=sqrt((i-xx).^2+(j-yy).^2+((k+0.5)-zz).^2);
            h=normpdf(d,0,sigma);
            J(1,1)=sum(sum(G11.*h));
            J(1,2)=sum(sum(G12.*h));
            J(1,3)=sum(sum(G13.*h));
            J(2,1)=sum(sum(G21.*h));
            J(2,2)=sum(sum(G22.*h));
            J(2,3)=sum(sum(G23.*h));
            J(3,1)=sum(sum(G31.*h));
            J(3,2)=sum(sum(G32.*h));
            J(3,3)=sum(sum(G33.*h));
            [V,D]=eig(J);
            Values(i-window,j-window,(k-1)*2+1,:,:)=D;
            Vectors(i-window,j-window,(k-1)*2+1,:,:)=V;
        end
    end
end

% V(:,1) 에 해당하는 eigen value D(1,1)
% V(:,2)==> D(2,2)
% weight1=

end