clear,clc;
format short
fid=fopen('E:\2000K\a-Fe-a-Fe_yure_lowquean_lowercooling\a-Fe-a-Fe_1800K\nvtchiyu.lammpstrj','rt'); 
for ii=1:1
c3=cell2mat(textscan(fid,'%f','headerlines',9));
b=0;
dd=zeros((length(c3)/8),8);
for i=1:8:length(c3)
    b=b+1;
  for  j=i:i+7
      k=j-i+1;
      dd(b,k)=c3(j,1);
  end
end
dd(dd(:,1)==0,:)=[];
d=zeros(size(dd,1),8);
for i=1:size(dd,1)
    for j=1:size(dd,1)
        if dd(j,1)==i
            d(i,:)=dd(j,:);
        end
    end
end
d(d(:,1)==0,:)=[];

ab_b=zeros(size(dd,1)*size(dd,1),18);
ab_bxyz=zeros(size(dd,1)*size(dd,1),8);
for i=1:size(d,1)
    cj=zeros(size(dd,1),8);
    for j=1:size(d,1)
        if d(i,1)~=d(j,1)&&abs(d(i,3)-d(j,3))<=3.325 && abs(d(i,4)-d(j,4))<=3.325 && abs(d(i,5)-d(j,5))<=3.325
            l=sqrt((d(i,3)-d(j,3))^2+(d(i,4)-d(j,4))^2+(d(i,5)-d(j,5))^2);
            if l<=3.325
                cj(j,:)=d(j,:);
            end
        end
    end
    cj(cj(:,1)==0,:)=[];
    if size(cj,1)>0
        num_ab(i,1)=i;
        num_ab(i,2)=size(cj,1);
        ab_b(i,1)=d(i,1);
        for k=1:size(cj,1)
            ab_b(i,k+1)=cj(k,1);
            ab_bxyz(sum(num_ab(:,2))-size(cj,1)+k,:)=cj(k,:);
        end
    end
    if size(cj,1)==0
        num_ab(i,1)=i;
        num_ab(i,2)=size(cj,1);
        ab_b(i,1)=d(i,1);
        ab_b(i,2)=0;
    end
end
ab_b(ab_b(:,1)==0,:)=[];
ab_bxyz(ab_bxyz(:,1)==0,:)=[];
B=zeros(10000,35);
ab_c=zeros(size(dd,1)*size(dd,1),18);
for i=1:size(ab_b,1)
    for j=2:size(ab_b,2)
        A=ab_b(i,:);
        if ab_b(i,j)>0
           B=ab_b(ab_b(i,j),:); 
        end
        C=zeros(1,size(ab_b,2));
        for k=2:size(A,2)
            for kk=2:size(B,2)
                if A(1,k)~=0 && B(1,kk)~=0 && A(1,k)==B(1,kk)
                    C(1,k)=A(1,k);
                end
            end
        end
        C(:,C(1,:)==0)=[];
        if size(C,2)>0
            ab_c((i-1)*18+j-1,1)=i;
            ab_c((i-1)*18+j-1,2)=ab_b(i,j);
            for kkk=1:size(C,2)
                ab_c((i-1)*18+j-1,2+kkk)=C(1,kkk);
            end
        end
    end
end
ab_c(ab_c(:,2)==0,:)=[];
for i=1:size(ab_c,1)
    t=0;
    tt=0;
    for j=3:size(ab_c,2)
        if ab_c(i,j)>0
            t=t+1;
        end
        for w=(j+1):size(ab_c,2)
            if ab_c(i,j)>0&&ab_c(i,w)>0
                l=sqrt((d(ab_c(i,j),3)-d(ab_c(i,w),3))^2+(d(ab_c(i,j),4)-d(ab_c(i,w),4))^2+(d(ab_c(i,j),5)-d(ab_c(i,w),5))^2);
                if l<=3.325
                    tt=tt+1;
                end
            end
        end
    end
    num_ab_c(i,1)=i;
    num_ab_c(i,2)=t;
    num_ab_c(i,3)=tt;
end
for i=1:6
    for j=1:6
        u=0;
        for k=1:size(num_ab_c,1)
            if num_ab_c(k,2)==i&&num_ab_c(k,3)==j
                u=u+1;
            end
        end
        num_jk((i-1)*6+j,1)=i;
        num_jk((i-1)*6+j,2)=j;
        num_jk((i-1)*6+j,3)=u;
    end
end
end
 fclose(fid);
 pinlv=sortrows(num_jk,3)