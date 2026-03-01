clear,clc;
format short
fid=fopen('E:\laser-cladding\c-Fe-c-Fe_yure\1\c-Fe-c-Fe_300K\4000down.lammpstrj','rt'); 
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
        if d(i,1)~=d(j,1)&&abs(d(i,3)-d(j,3))<=2.9 && abs(d(i,4)-d(j,4))<=2.9 && abs(d(i,5)-d(j,5))<=2.9  
            l=sqrt((d(i,3)-d(j,3))^2+(d(i,4)-d(j,4))^2+(d(i,5)-d(j,5))^2);
            if l<=2.9
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
lll=ab_b;
%AA=zeros(size(ab_b,1)*(size(ab_b,2)-1),8);
lll(:,1)=[];
bb=lll';
AA=bb(:);
AA(AA(:,1)==0,:)=[]
CC=zeros(size(ab_b,1)*(size(ab_b,2)),8);
end
for qq=1:9216    
        for uu=1:17
          if lll(qq,uu)>0
            tt=(qq*17-16)+uu;
            CC(tt,:)=d(qq,:);
          end
        end
end
CC(CC(:,1)==0,:)=[];
length=zeros(100220,1);
for gg=1:100220
    hh=gg;
        length(gg,1)=sqrt((CC(gg,3)-ab_bxyz(hh,3))^2+(CC(gg,4)-ab_bxyz(hh,4))^2+(CC(gg,5)-ab_bxyz(hh,5))^2);
end         

 