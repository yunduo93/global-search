function f=mean_negtive(x)
[h,w]=size(x);
num=h*w;
i=(num+1)/2;
s=0;
ind=0;
if(x(i)<0)
    for j=1:num
     if(x(j)>0)
       s=s+x(j);
       ind=ind+1;
     end
    end
    if(ind~=0)
        f=s/ind;
    end
   
    if(s==0||ind==0)
        f=0.01;
    end
    
else
    f=x(i);
end
