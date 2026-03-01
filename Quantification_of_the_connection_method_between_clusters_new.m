
clear,clc;
format short
fid=fopen('E:\laser-cladding\c-Fe-c-Fe_coolingrate\c-Fe-c-Fe_300K\nvtchiyu.lammpstrj','rt'); 
jl=3.325;
q1=101;
tuancu=zeros(9216,20);
for ii=1:q1
    c3=cell2mat(textscan(fid,'%f','headerlines',9));
    if ii==101
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
        tuancu(:,1)=d(:,1);
        for k=1:size(d,1)
            number=0;
            for i=1:size(d,1)
                if d(i,1)~=d(k,1)&&abs(d(i,3)-d(k,3))<=jl && abs(d(i,4)-d(k,4))<=jl && abs(d(i,5)-d(k,5))<=jl&&d(i,3)~=0&&d(k,3)~=0&&d(i,4)~=0&&d(k,4)~=0&&d(i,5)~=0&&d(k,5)~=0  
                    l=sqrt((d(i,3)-d(k,3))^2+(d(i,4)-d(k,4))^2+(d(i,5)-d(k,5))^2);
                    if l<=jl
                        number=number+1;
                         tuancu(k,number+1)=d(i,1);
                    end
                end
            end
        end
    end
end
guanlian_1=zeros(9216,20);guanlian_2=zeros(9216,20);guanlian_3=zeros(9216,20);guanlian_4=zeros(9216,20);guanlian_5=zeros(9216,20);
g1=0;g2=0;g3=0;g4=0;g5=0;
for i=1:size(tuancu,1)
    n1=0;n2=0;n3=0;n4=0;n5=0;
    for j=1:size(tuancu,1)
        tt=0;ttt=0;
       for k=2:size(tuancu,2)
            for p=2:size(tuancu,2)
                if tuancu(i,k)==tuancu(j,p)&&tuancu(i,k)~=0&&tuancu(j,p)~=0
                    if (sqrt((d(i,3)-d(j,3))^2+(d(i,4)-d(j,4))^2+(d(i,5)-d(j,5))^2)>jl)&&(sqrt((d(i,3)-d(j,3))^2+(d(i,4)-d(j,4))^2+(d(i,5)-d(j,5))^2)<3*jl)
                        tt=tt+1;
                    else
                        ttt=ttt+1;
                    end
                end
            end
       end
       switch tt
           case 1
                n1=n1+1;g1=g1+1;
                guanlian_1(i,1)=tuancu(i,1);guanlian_1(i,n1+1)=tuancu(j,1); 
            case 2
                n2=n2+1;g2=g2+1;
                guanlian_2(i,1)=tuancu(i,1);guanlian_2(i,n2+1)=tuancu(j,1); 
            case 3
                n3=n3+1;g3=g3+1;
                guanlian_3(i,1)=tuancu(i,1);guanlian_3(i,n3+1)=tuancu(j,1);
            case {4,5,6,7,8,9}
                n4=n4+1;g4=g4+1;
                guanlian_4(i,1)=tuancu(i,1);guanlian_4(i,n4+1)=tuancu(j,1);
       end
        if ttt~=0
            n5=n5+1;g5=g5+1;
            guanlian_5(i,1)=tuancu(i,1);guanlian_5(i,n5+1)=tuancu(j,1);
        end
    end
end
guanlian_1(guanlian_1(:,1)==0,:)=[];guanlian_2(guanlian_2(:,1)==0,:)=[];guanlian_3(guanlian_3(:,1)==0,:)=[];guanlian_4(guanlian_4(:,1)==0,:)=[];guanlian_5(guanlian_5(:,1)==0,:)=[];
tongji=zeros(5,2);
tongji(1,1)=1;tongji(2,1)=2;tongji(3,1)=3;tongji(4,1)=4;tongji(5,1)=5;
tongji(1,2)=g1/2;tongji(2,2)=g2/2;tongji(3,2)=g3/2;tongji(4,2)=g4/2;tongji(5,2)=g5/2;
fclose(fid);