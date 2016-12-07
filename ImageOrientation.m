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
for k=1:size(image,3)
    [G1,G2]=imgradientxy(PadImage(:,:,k));
    if k==size(image,3)
        G3=PadImage(:,:,k)-PadImage(:,:,1);
    else
        G3=PadImage(:,:,k)-PadImage(:,:,k+1);
    end
    for i=window+1:size(PadImage,1)-window
        for j=window+1:size(PadImage,2)-window
            ii=i-window+1:i+window;
            jj=j-window+1:j+window;
            [xx,yy]=meshgrid(ii,jj);
            d=sqrt((i-xx).^2+(j-yy).^2);
            J(1,1)=sum(sum((G1(ii,jj).*G1(ii,jj)).*normpdf(d,0,sigma)));
            J(1,2)=sum(sum((G1(ii,jj).*G2(ii,jj)).*normpdf(d,0,sigma)));
            J(1,3)=sum(sum((G1(ii,jj).*G3(ii,jj)).*normpdf(d,0,sigma)));
            J(2,1)=sum(sum((G2(ii,jj).*G1(ii,jj)).*normpdf(d,0,sigma)));
            J(2,2)=sum(sum((G2(ii,jj).*G2(ii,jj)).*normpdf(d,0,sigma)));
            J(2,3)=sum(sum((G2(ii,jj).*G3(ii,jj)).*normpdf(d,0,sigma)));
            J(3,1)=sum(sum((G3(ii,jj).*G1(ii,jj)).*normpdf(d,0,sigma)));
            J(3,2)=sum(sum((G3(ii,jj).*G2(ii,jj)).*normpdf(d,0,sigma)));
            J(3,3)=sum(sum((G3(ii,jj).*G3(ii,jj)).*normpdf(d,0,sigma)));
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