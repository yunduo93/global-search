clc;
clear;
addpath('./BM3D');
Path='./images/';
num_name=2;            %1, 2, 3
str_name_i='_I.tiff';  %.tiff, .bmp
str_name_o='_O.tiff';
I=imread([Path,num2str(num_name),str_name_i]);   
O=imread([Path,num2str(num_name),str_name_o]);
infc=load(['./',num2str(num_name),'.mat']);
I=im2double(I);
O=im2double(O);
[col,row,n]=size(I);
num=col*row;

Re_min=O-I;
Re_max=O+I;

infty=solve_A(Re_min,Re_max)*infc.mn(2);

lb = 0.2+infc.mn(4);
ub = 0.5;
x0 = [0.2 0.3 0.6 0.8 1];
for i=1:3
    I1=im2double(im2uint8((Re_min(:,:,i))));
    I2=im2double(im2uint8((Re_max(:,:,i))));
    I1=I1(:);
    I2=I2(:);
    a=[I1,I2];
    a=sortrows(a,1);
    a=round(a*100)/100;
    j=1;
    while(j<255*255)
       if(a(j,1)~=0)
          break;
       end
       j=j+1;
    end
    a([1:j],:)=[];
    [b m n]=unique(a,'rows');
    c=tabulate(n');
    s=[a(m(c(:,1)),:),c(:,2:end)];
    [h,w]=size(s);
    
    syms p f
    k=Re_min(:,:,i);
    m=Re_max(:,:,i);
    E_k=sum(k(:))/num;
    f0=0;
    f1=0;
    for o=1:h
        k=s(o,1);
        m=s(o,2);
        coeff=s(o,4);
        f0=(k.*(p*m-k)./(p*infty(i)-k))*coeff;
        f1=E_k*(p*m-k)./(p*infty(i)-k)*coeff+f1-f0;
    end
    f=abs(f1);
    G=inline(f,'p');
    %%
     problem = createOptimProblem('fmincon','objective',G,'x0',x0(1),'lb',lb,...
            'ub',ub,'options',optimset('Algorithm','SQP','Disp','none'));
    gs = GlobalSearch;
    xgs = run(gs,problem);
    %%
     P(i)=xgs;
     ub=xgs;
end

p_max=0;
for i=1:3
    if(P(i)>p_max)
        p_max=P(i);
    end
end
for i=1:3
    if(abs(P(i)-p_max)>0.05)
      P(i)=p_max-0.018;
    end
end 

P=P*infc.mn(3);
e=1.3;
e1=1;

 %% 
 A=zeros(col,row,3);
 w=zeros(col,row,3);
 L=zeros(col,row,3);
 D_test=zeros(col,row,3);
 L1=zeros(col,row,3);
 L2=zeros(col,row,3);
 trans=zeros(col,row,3);
  
 for i=1:3
    A(:,:,i)=Re_min(:,:,i)./P(i);
 end
  D_test=Re_max-A;  

%% 
 for i=1:3 
     TEMP=A(:,:,i)./(infty(i)*e1); 
     t=1-TEMP;
      while(min(t(:))<0)
        t=nlfilter(t,[5 5],@mean_negtive);
      end
    [PSNR,t] = BM3D(1, t, 10 , 'lc', 1);
     t=min(t,1);
     trans(:,:,i)=t;
     D=D_test(:,:,i);
     iterator=0;
     while(min(D(:))<0)
        D=nlfilter(D,[5 5],@mean_negtive);
        iterator=iterator+1;
        if(iterator==5)
            disp('Over max iterator!');
            break;
        end
     end
    
     L1(:,:,i)=D./t; 
     L1(:,:,i)=im2double(im2uint16(L1(:,:,i)));
 end
 
 %% 
P=P*e;
 for i=1:3 
      A(:,:,i)=Re_min(:,:,i)./P(i);
     TEMP=A(:,:,i)./(infty(i)*e1); 
     t=1-TEMP;
      while(min(t(:))<0)
        t=nlfilter(t,[5 5],@mean_negtive);
      end
     [PSNR,t] = BM3D(1, t, 10 , 'lc', 1);
     t=min(t,1);
    
    D=Re_max(:,:,i)-A(:,:,i);  
     iterator=0;
     while(min(D(:))<0)
        D=nlfilter(D,[5 5],@mean_negtive);
        iterator=iterator+1;
        if(iterator==5)
            disp('Over max iterator!');
            break;
        end
     end
     L2(:,:,i)=D./t;
     L2(:,:,i)=im2double(im2uint16(L2(:,:,i)));
 end
 
 for i=1:3
     [L2(:,:,i)] = BM3D_CFA(L2(:,:,i), 25);
 end
 
 H=L1-L2;
 temp=H.^2;
 for i=1:3
    temp_bw=im2bw(imadjust(temp(:,:,i)),0.15);
    w(:,:,i)=1-trans(:,:,1);
    w(:,:,i)=imadjust(w(:,:,i));
    H(:,:,i)=H(:,:,i).*temp_bw.*w(:,:,i);
 end
 
 K=L1-H;
 for i=1:3
    [K(:,:,i)] = BM3D_CFA(K(:,:,i), 15);
    L(:,:,i)=L1(:,:,i).*(1-w(:,:,i))+K(:,:,i).*w(:,:,i);
 end
 
 if(infc.mn(1)==1)
 gray=rgb2gray(L);
 seg=im2bw(gray,0.6);

cel=[1.6,1.4,1.3];
bw=ones(col,row)*1.6;

se=strel('square',5);
bw1=imdilate(seg,se);
bw2=imdilate(bw1,se);
bw3=imdilate(bw2,se);

for i=1:3
seg1=(1-bw3)*cel(i);
seg2=(bw3-bw2)*(cel(i)*0.9);
seg3=(bw2-bw1)*(cel(i)*0.8);
seg4=bw1*(cel(i)*0.7);

mag=seg1+seg2+seg3+seg4;

L(:,:,i)=L(:,:,i).*mag;
end

seg1=(1-bw3)*1.3;
seg2=(bw3-bw2)*1.5;
seg3=(bw2-bw1)*1.3;
seg4=bw1*1.1;

mag=seg1+seg2+seg3+seg4;
for i=1:3
    L(:,:,i)=L(:,:,i).*mag;
end
 end
imshow(L);

