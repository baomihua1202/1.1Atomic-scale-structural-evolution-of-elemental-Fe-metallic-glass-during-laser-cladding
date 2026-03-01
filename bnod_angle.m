clear,clc;
format short
fid=fopen('E:\2000K\a-Fe-a-Fe_yure_lowquean_highercooling\a-Fe-a-Fe_1800K\nvtchiyu.lammpstrj','rt'); 
for ii=1:1
c3=cell2mat(textscan(fid,'%f','headerlines',9));
b=0;
dd=zeros((length(c3)/8),8);
for i=1:8:length(c3)
    b=b+1;
  for  j=i:i+7
      k=j-i+1;
      dd(b,k)=c3(j,1);%dd中数据混乱排布 
  end
end
dd(dd(:,1)==0,:)=[];
d=zeros(size(dd,1),8);
for i=1:size(dd,1)
    for j=1:size(dd,1)
        if dd(j,1)==i
            d(i,:)=dd(j,:);%将数据进行排序
        end
    end
end
d(d(:,1)==0,:)=[];

%第一步：确定A原子，找与A原子成A-B键的B原子ID、坐标信息等
%A原子ID为1-9216
ab_b=zeros(10000,18);%找成A-B键的B原子ID，所有原子成键的ID都找到了,第一列为A原子ID，2-18列为B原子ID
ab_bxyz=zeros(200000,8);%将找到的成A-B键的B原子信息进行存储
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
end
%kkk=ab_b;
%for qq=1:9216
   % for ww=2:16
       % if kkk(qq,ww)==0
          %  ab_b(qq,ww)=i;
       % end
    %end
%end
ab_b(ab_b(:,1)==0,:)=[];
ab_bxyz(ab_bxyz(:,1)==0,:)=[];
K=zeros(size(ab_b,1),300);
for z=1:9216
    for x=2:16
        for c=(x+1):17
          if ab_b(z,x)>0&&ab_b(z,c)>0&&ab_b(ab_b(z,x),c)~=d(z,1);     %ab_b(i,j)~=ab_b(ab_b(i,j),k)&&(ab_b(ab_b(i,j),k))>0&&ab_b(ab_b(i,j),k)~=d(i,1)
            a=(d(ab_b(z,1),3)-d(ab_b(z,x),3))^2+(d(ab_b(z,1),4)-d(ab_b(z,x),4))^2+(d(ab_b(z,1),5)-d(ab_b(z,x),5))^2;
            b=(d(ab_b(z,1),3)-d(ab_b(z,c),3))^2+(d(ab_b(z,1),4)-d(ab_b(z,c),4))^2+(d(ab_b(z,1),5)-d(ab_b(z,c),5))^2;
            a1=sqrt(a);
            b1=sqrt(b);
            ab=(d(ab_b(z,1),3)-d(ab_b(z,x),3)).*(d(ab_b(z,1),3)-d(ab_b(z,c),3))+(d(ab_b(z,1),4)-d(ab_b(z,x),4)).*(d(ab_b(z,1),4)-d(ab_b(z,c),4))+(d(ab_b(z,1),5)-d(ab_b(z,x),5)).*(d(ab_b(z,1),5)-d(ab_b(z,c),5));
            K(z,(17*x-33+c))=(acos(ab/(a1*b1)))*180/pi;
          end
        end
    end
end
K(find(K==0))=[];
L=K';
M=round(L);
N=tabulate(M)