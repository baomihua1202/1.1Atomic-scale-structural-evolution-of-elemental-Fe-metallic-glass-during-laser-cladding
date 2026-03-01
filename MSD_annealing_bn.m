clear,clc;
format short
fid4=fopen('E:\2000K\a-Fe-a-Fe_yure_lowquean_highercooling\a-Fe-a-Fe_650K\2000down.lammpstrj','rt'); 
fid1=fopen('E:\2000K\a-Fe-a-Fe_yure_lowquean_highercooling\a-Fe-a-Fe_650K\annealing1.lammpstrj','rt');
fid2=fopen('E:\2000K\a-Fe-a-Fe_yure_lowquean_highercooling\a-Fe-a-Fe_650K\650chiyu1.lammpstrj','rt'); 
fid3=fopen('E:\2000K\a-Fe-a-Fe_yure_lowquean_highercooling\a-Fe-a-Fe_650K\annealing2.lammpstrj','rt'); 
fid11=fopen('E:\2000K\a-Fe-a-Fe_yure_lowquean_highercooling\a-Fe-a-Fe_650K\annealing.txt','w');
wd=650;q2=101;
for i=1:1
    if wd==650
        q1=55;q3=15;
    end
    if wd==1000
        q1=41;q3=29;
    end
    if wd==1400
       q1=25;q3=45;
    end
    if wd==1800
        q1=9;q3=61;
    end
end
for ii=1:1
    c3=cell2mat(textscan(fid4,'%f','headerlines',9));
    b=0;
    dd0=zeros((length(c3)/3),3);
    for i=1:8:length(c3)
        b=b+1;
      for  j=i:i+7
          k=j-i+1;
          dd0(b,k)=c3(j,1);
      end
    end
    dd0(dd0(:,1)==0,:)=[]; 
    d0=zeros(size(dd0,1),8);
    for i=1:size(dd0,1)
        for j=1:size(dd0,1)
            if dd0(j,1)==i
                d0(i,:)=dd0(j,:);
            end
        end
    end
    d0(d0(:,1)==0,:)=[]; 
end

msd2_1=zeros(q1+1,9216+9216);msd3_1=zeros(q1+1,9216+9216);msd4_1=zeros(q1+1,9216);
for ii=1:q1
    c3=cell2mat(textscan(fid1,'%f','headerlines',9));
    b=0;
    dd1=zeros((length(c3)/3),3);
    for i=1:8:length(c3)
        b=b+1;
      for  j=i:i+7
          k=j-i+1;
          dd1(b,k)=c3(j,1);
      end
    end
    dd1(dd1(:,1)==0,:)=[]; 
    d1=zeros(size(dd1,1),8);
    for i=1:size(dd1,1)
        for j=1:size(dd1,1)
            if dd1(j,1)==i
                d1(i,:)=dd1(j,:);
            end
        end
    end
    d1(d1(:,1)==0,:)=[];
    for i=1:size(d1,1)
        if d1(i,2)==2
            l2=(d1(i,3)-d0(i,3))^2+(d1(i,4)-d0(i,4))^2+(d1(i,5)-d0(i,5))^2;
            msd2_1(1,i)=d1(i,1);
            msd2_1(ii+1,i)=l2;
        end
        if d1(i,2)==3
            l23=(d1(i,3)-d0(i,3))^2+(d1(i,4)-d0(i,4))^2+(d1(i,5)-d0(i,5))^2;
            msd3_1(1,i)=d1(i,1);
            msd3_1(ii+1,i)=l23;
        end
        if d1(i,2)==4
            l4=(d1(i,3)-d0(i,3))^2+(d1(i,4)-d0(i,4))^2+(d1(i,5)-d0(i,5))^2;
            msd4_1(1,i)=d1(i,1);
            msd4_1(ii+1,i)=l4;
        end
    end
end
msd2_1(:,msd2_1(1,:)==0)=[];msd3_1(:,msd3_1(1,:)==0)=[]; msd4_1(:,msd4_1(1,:)==0)=[]; 
msd234_1=zeros(q1,size(msd3_1,2)+size(msd4_1,2)+size(msd2_1,2));
msd34_1=zeros(q1,size(msd3_1,2)+size(msd4_1,2));
msd234_1=[msd2_1,msd3_1,msd4_1];
msd34_1=[msd3_1,msd4_1];
a234_1=mean(msd234_1,2);a34_1=mean(msd34_1,2);a3_1=mean(msd3_1,2);a4_1=mean(msd4_1,2);
a234_1(1,:)=[];a34_1(1,:)=[];a3_1(1,:)=[];a4_1(1,:)=[];
a_1=zeros(q1,5);
a_1(:,1)=linspace(1,q1,q1);a_1(:,2)=a234_1;a_1(:,3)=a34_1;a_1(:,4)=a3_1;a_1(:,5)=a4_1;

msd2_2=zeros(q2+1,9216+9216);msd3_2=zeros(q2+1,9216+9216);msd4_2=zeros(q2+1,9216);
for ii=1:q2
    c3=cell2mat(textscan(fid2,'%f','headerlines',9));
    b=0;
    dd2=zeros((length(c3)/3),3);
    for i=1:8:length(c3)
        b=b+1;
      for  j=i:i+7
          k=j-i+1;
          dd2(b,k)=c3(j,1);
      end
    end
    dd2(dd2(:,1)==0,:)=[]; 
    d2=zeros(size(dd2,1),8);
    for i=1:size(dd2,1)
        for j=1:size(dd2,1)
            if dd2(j,1)==i
                d2(i,:)=dd2(j,:);
            end
        end
    end
    d2(d2(:,1)==0,:)=[];
    for i=1:size(d2,1)
        if d2(i,2)==2
            l2=(d2(i,3)-d0(i,3))^2+(d2(i,4)-d0(i,4))^2+(d2(i,5)-d0(i,5))^2;
            msd2_2(1,i)=d2(i,1);
            msd2_2(ii+1,i)=l2;
        end
        if d2(i,2)==3
            l23=(d2(i,3)-d0(i,3))^2+(d2(i,4)-d0(i,4))^2+(d2(i,5)-d0(i,5))^2;
            msd3_2(1,i)=d2(i,1);
            msd3_2(ii+1,i)=l23;
        end
        if d2(i,2)==4
            l4=(d2(i,3)-d0(i,3))^2+(d2(i,4)-d0(i,4))^2+(d2(i,5)-d0(i,5))^2;
            msd4_2(1,i)=d2(i,1);
            msd4_2(ii+1,i)=l4;
        end
    end
end
msd2_2(:,msd2_2(1,:)==0)=[];msd3_2(:,msd3_2(1,:)==0)=[]; msd4_2(:,msd4_2(1,:)==0)=[]; 
msd234_2=zeros(q2,size(msd3_2,2)+size(msd4_2,2)+size(msd2_2,2));
msd34_2=zeros(q2,size(msd3_2,2)+size(msd4_2,2));
msd234_2=[msd2_2,msd3_2,msd4_2];
msd34_2=[msd3_2,msd4_2];
a234_2=mean(msd234_2,2);a34_2=mean(msd34_2,2);a3_2=mean(msd3_2,2);a4_2=mean(msd4_2,2);
a234_2(1,:)=[];a34_2(1,:)=[];a3_2(1,:)=[];a4_2(1,:)=[];
a_2=zeros(q2,5);
a_2(:,1)=linspace(1,q2,q2);a_2(:,2)=a234_2;a_2(:,3)=a34_2;a_2(:,4)=a3_2;a_2(:,5)=a4_2;
msd2_3=zeros(q3+1,9216+9216);msd3_3=zeros(q3+1,9216+9216);msd4_3=zeros(q3+1,9216);
for ii=1:q3
    c3=cell2mat(textscan(fid3,'%f','headerlines',9));
    b=0;
    dd3=zeros((length(c3)/3),3);
    for i=1:8:length(c3)
        b=b+1;
      for  j=i:i+7
          k=j-i+1;
          dd3(b,k)=c3(j,1);
      end
    end
    dd3(dd3(:,1)==0,:)=[]; 
    d3=zeros(size(dd3,1),8);
    for i=1:size(dd3,1)
        for j=1:size(dd3,1)
            if dd3(j,1)==i
                d3(i,:)=dd3(j,:);
            end
        end
    end
    d3(d3(:,1)==0,:)=[];
    for i=1:size(d3,1)
        if d3(i,2)==2
            l2=(d3(i,3)-d0(i,3))^2+(d3(i,4)-d0(i,4))^2+(d3(i,5)-d0(i,5))^2;
            msd2_3(1,i)=d3(i,1);
            msd2_3(ii+1,i)=l2;
        end
        if d3(i,2)==3
            l23=(d3(i,3)-d0(i,3))^2+(d3(i,4)-d0(i,4))^2+(d3(i,5)-d0(i,5))^2;
            msd3_3(1,i)=d3(i,1);
            msd3_3(ii+1,i)=l23;
        end
        if d3(i,2)==4
            l4=(d3(i,3)-d0(i,3))^2+(d3(i,4)-d0(i,4))^2+(d3(i,5)-d0(i,5))^2;
            msd4_3(1,i)=d3(i,1);
            msd4_3(ii+1,i)=l4;
        end
    end
end
msd2_3(:,msd2_3(1,:)==0)=[];msd3_3(:,msd3_3(1,:)==0)=[]; msd4_3(:,msd4_3(1,:)==0)=[]; 
msd234_3=zeros(q3,size(msd3_3,2)+size(msd4_3,2)+size(msd2_3,2));
msd34_3=zeros(q3,size(msd3_3,2)+size(msd4_3,2));
msd234_3=[msd2_3,msd3_3,msd4_3];
msd34_3=[msd3_3,msd4_3];
a234_3=mean(msd234_3,2);a34_3=mean(msd34_3,2);a3_3=mean(msd3_3,2);a4_3=mean(msd4_3,2);
a234_3(1,:)=[];a34_3(1,:)=[];a3_3(1,:)=[];a4_3(1,:)=[];
a_3=zeros(q3,5);
a_3(:,1)=linspace(1,q3,q3);a_3(:,2)=a234_3;a_3(:,3)=a34_3;a_3(:,4)=a3_3;a_3(:,5)=a4_3;
a=zeros(q1+q2+q3,5);
a=[a_1;a_2;a_3];
a(:,1)=linspace(1,size(a,1),size(a,1));
for i=1:size(a,1)
    fprintf(fid11,'%d %f %f %f %f\n',i,a(i,2),a(i,3),a(i,4),a(i,5));
end