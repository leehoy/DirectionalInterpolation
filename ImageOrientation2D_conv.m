function [ J11,J12,J21,J22] = ImageOrientation2D_conv( image,sigma,window)
%Image orientation function for 2d sinogram case
% z-direction in omitted
% [G1,G2]=imgradientxy(image(:,:,k));
% G3=image(:,:,k)-image(:,:,k+1);
% image=image';
% Vectors=zeros(size(image,2),size(image,1),2,2);
% Values=zeros(size(image,2),size(image,1),2,2);
% PadImage=zeros(size(image,1)+2*window,size(image,2)+2*window);
% PadImage=padarray(image,[window window],'circular','both');
[G1,G2]=imgradientxy(image); % G1 - x direction, G2 - y direction
% x=-window+1:size(image,2)+window;
% y=-window+1:size(image,1)+window;
x=-window:window;
y=-window:window;
[xx,yy]=meshgrid(x,y);
G11=G1.*G1;
G12=G1.*G2;
G21=G2.*G1;
G22=G2.*G2;
d=sqrt((0-xx).^2+(0.5-yy).^2);
h=normpdf(d,0,sigma);
J11=conv2(G11,h,'same');
J12=conv2(G12,h,'same');
J21=conv2(G21,h,'same');
J22=conv2(G22,h,'same');
% for i=1:size(image,1)
%     for j=2:size(image,2)*2
% for i=window+1:size(PadImage,1)-window
%     for j=window+2:2:size(PadImage,2)*2-window
%         ii=i-window+1:i+window;
%         jj=j-window+1:j+window;
%         [xx,yy]=meshgrid(ii,jj);
%         d=sqrt((i-xx).^2+(j-yy).^2);
%         J(1,1)=sum(sum((G1(ii,jj).*G1(ii,jj)).*normpdf(d,0,sigma)));
%         J(1,2)=sum(sum((G1(ii,jj).*G2(ii,jj)).*normpdf(d,0,sigma)));
% %             J(1,3)=sum(sum((G1(ii,jj).*G3(ii,jj)).*normpdf(d,0,sigma)));
%         J(2,1)=sum(sum((G2(ii,jj).*G1(ii,jj)).*normpdf(d,0,sigma)));
%         J(2,2)=sum(sum((G2(ii,jj).*G2(ii,jj)).*normpdf(d,0,sigma)));
% %             J(2,3)=sum(sum((G2(ii,jj).*G3(ii,jj)).*normpdf(d,0,sigma)));
% %             J(3,1)=sum(sum((G3(ii,jj).*G1(ii,jj)).*normpdf(d,0,sigma)));
% %             J(3,2)=sum(sum((G3(ii,jj).*G2(ii,jj)).*normpdf(d,0,sigma)));
% %             J(3,3)=sum(sum((G3(ii,jj).*G3(ii,jj)).*normpdf(d,0,sigma)));
%         [V,D]=eig(J);
%         Values(i-window,j-window,:,:)=D;
%         Vectors(i-window,j-window,:,:)=V;
%     end
% end
% % % [G1,G2]=imgradientxy(PadImage);
% % 
% % %     if k==size(image,3)
% % %         G3=PadImage(:,:,k)-PadImage(:,:,1);
% % %     else
% % %         G3=PadImage(:,:,k)-PadImage(:,:,k+1);
% % %     end
% % for i=window+1:size(PadImage,1)-window
% %     for j=window+1:size(PadImage,2)-window
% %         ii=i-window+1:i+window;
% %         jj=j-window+1:j+window;
% %         [xx,yy]=meshgrid(ii,jj);
% %         d=sqrt((i-xx).^2+(j-yy).^2);
% %         J(1,1)=sum(sum((G1(ii,jj).*G1(ii,jj)).*normpdf(d,0,sigma)));
% %         J(1,2)=sum(sum((G1(ii,jj).*G2(ii,jj)).*normpdf(d,0,sigma)));
% % %             J(1,3)=sum(sum((G1(ii,jj).*G3(ii,jj)).*normpdf(d,0,sigma)));
% %         J(2,1)=sum(sum((G2(ii,jj).*G1(ii,jj)).*normpdf(d,0,sigma)));
% %         J(2,2)=sum(sum((G2(ii,jj).*G2(ii,jj)).*normpdf(d,0,sigma)));
% % %             J(2,3)=sum(sum((G2(ii,jj).*G3(ii,jj)).*normpdf(d,0,sigma)));
% % %             J(3,1)=sum(sum((G3(ii,jj).*G1(ii,jj)).*normpdf(d,0,sigma)));
% % %             J(3,2)=sum(sum((G3(ii,jj).*G2(ii,jj)).*normpdf(d,0,sigma)));
% % %             J(3,3)=sum(sum((G3(ii,jj).*G3(ii,jj)).*normpdf(d,0,sigma)));
% %         [V,D]=eig(J);
% %         Values(i-window,j-window,:,:)=D;
% %         Vectors(i-window,j-window,:,:)=V;
% %     end
% % end
% 
% % V(:,1) 에 해당하는 eigen value D(1,1)
% % V(:,2)==> D(2,2)
% % weight1=

end