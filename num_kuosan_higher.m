clear,clc;
format short
fid=fopen('F:\laser-cladding\c-Fe-c-Fe_yure\1\c-Fe-c-Fe_1800K\4000down.lammpstrj','rt');
for ii=1:1
c3=cell2mat(textscan(fid,'%f','headerlines',9));
b=0;
d=zeros((length(c3)/8),8);
for i=1:8:length(c3)
    b=b+1;
  for  j=i:i+7
      k=j-i+1;
      d(b,k)=c3(j,1);
  end
end
d(d(:,1)==0,:)=[];
jidi=zeros(size(d,1),8);
for i=1:size(d,1)
    if d(i,2)==2||d(i,2)==3
        jidi(i,:)=d(i,:);
    end
end
jidi(jidi(:,1)==0,:)=[];
z_max=max(jidi(:,5));
z_max=z_max+0.2;
end
 fclose(fid);
fid1=fopen('F:\laser-cladding\c-Fe-c-Fe_coolingrate\c-Fe-c-Fe_1800K\4000chiyu1.lammpstrj','rt'); 
fid2=fopen('F:\laser-cladding\c-Fe-c-Fe_coolingrate\c-Fe-c-Fe_1800K\n_ks_4000chiyu1-2.data','w');
fid3=fopen('F:\laser-cladding\c-Fe-c-Fe_coolingrate\c-Fe-c-Fe_1800K\Fe_ks_x_4000chiyu1-2.xyz','w');
fid4=fopen('F:\laser-cladding\c-Fe-c-Fe_coolingrate\c-Fe-c-Fe_1800K\Fe_ks_s_4000chiyu1-2.xyz','w');
 for jj=1:401 
c3=cell2mat(textscan(fid1,'%f','headerlines',9));
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

jidi1=zeros(size(dd,1),8);
rfc=zeros(size(dd,1),8);
for i=1:size(dd,1)
    if dd(i,2)==2||dd(i,2)==3
        jidi1(i,:)=dd(i,:);
    end
    if dd(i,2)==4
        rfc(i,:)=dd(i,:);
    end
end
jidi1(jidi1(:,1)==0,:)=[];
rfc(rfc(:,1)==0,:)=[];
ks_x=zeros(size(dd,1),8);
for i=1:size(jidi1,1)
    if jidi1(i,5)>z_max
        ks_x(i,:)=jidi1(i,:);
    end
end
ks_x(ks_x(:,1)==0,:)=[];

ks_s=zeros(size(dd,1),8);
for i=1:size(rfc,1)
    if rfc(i,5)<z_max
        ks_s(i,:)=rfc(i,:);
    end
end
ks_s(ks_s(:,1)==0,:)=[];

num_ks(jj,1)=jj;
num_ks(jj,2)=size(ks_x,1);
num_ks(jj,3)=size(ks_s,1);

fprintf(fid2,'%d %d %d\n',num_ks(jj,1),num_ks(jj,2),num_ks(jj,3));

if num_ks(jj,2)<1
    cc=1;
    fprintf(fid3,'%d\n',cc);  
    c2='Atoms. Timestep: nn';
    fprintf(fid3,'%s\n',c2); 
    fprintf(fid3,'%s %d %d %d\n','Fe',18.5,18.5,18.5);
end

if num_ks(jj,2)>0
    cc=num_ks(jj,2);
    fprintf(fid3,'%d\n',cc);  
    c2='Atoms. Timestep: nn';
    fprintf(fid3,'%s\n',c2); 
    for i=1:size(ks_x,1)
        fprintf(fid3,'%s %d %d %d\n','Fe',ks_x(i,3),ks_x(i,4),ks_x(i,5));
    end
end

if num_ks(jj,3)<1
    cc=1;
    fprintf(fid4,'%d\n',cc);  
    c2='Atoms. Timestep: nn';
    fprintf(fid4,'%s\n',c2); 
    fprintf(fid4,'%s %d %d %d\n','Fe',38.5,38.5,38.5);
end

if num_ks(jj,3)>0
    cc=num_ks(jj,3);
    fprintf(fid4,'%d\n',cc);  
    c2='Atoms. Timestep: nn';
    fprintf(fid4,'%s\n',c2); 
    for i=1:size(ks_s,1)
        fprintf(fid4,'%s %d %d %d\n','Fe',ks_s(i,3),ks_s(i,4),ks_s(i,5));
    end
end
 end
   fclose(fid1);
   fclose(fid2);
   fclose(fid3);
   fclose(fid4);
   
   