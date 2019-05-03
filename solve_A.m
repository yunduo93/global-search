function air_light = solve_A(I_min,I_total)
 
[m,n,~]=size(I_total);
imsize=m*n;
bw=I_min;
numpx=floor(imsize*0.05);
J=reshape(I_min,imsize,3);
a=zeros(1,3);

for i=1:3
    [J(:,i),indices]=sort(J(:,i));
    T=J(:,i);
    a(i)=T(indices(round(0.001*numpx)));
    I_minc=bw(:,:,i);
    I_minc(I_minc>a(i))=1;
    I_minc(I_minc<=a(i))=0;
    bw(:,:,i)=I_minc;
end

for i=1:3
    I_res(:,:,i)=I_total(:,:,i).*bw(:,:,i);
end
J=reshape(I_res,imsize,3);
temp=I_res;

for i=1:3
   [J(:,i),indices]=sort(J(:,i));
    T=J(:,i);
    a(i)=T(indices(round(0.001*numpx)));
    temp(temp>a(i))=1;
    temp(temp<=a(i))=0;
    I_res(:,:,i)=I_res(:,:,i).*temp(:,:,i);
    S=sum(sum(I_res(:,:,i)~=0,2));
    air_light(i)=sum(sum(I_res(:,:,i)))/S;
end


% for i=1:3
%     I_result(:,:,i)=bw(:,:,i).*I_total(:,:,i);
%     test=I_result(:,:,i);
%     S=sum(sum(I_result(:,:,i)~=0,2));
%     air_light(i)=sum(sum(I_result(:,:,i)))/S;
% end


