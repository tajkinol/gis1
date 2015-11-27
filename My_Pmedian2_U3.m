clear all;clc;
%[num,txt,raw] = xlsread('Matr_1.xls','Лист1');
[num,txt,raw] = xlsread('Улуг-Хем.xls','Лист1');

M=num;% матрица расстояний
n=size(M,1);%число узлов
nodes=1:n; %номера узлов
p=2;%число медиан
c=nchoosek(nodes,p);% все сочетания
R_c=c;
R_c(:,size(c,2)+1)=0;
chislo_soch=size(c,1);
%nc_vse=zeros(chislo_soch,n);
for kk=1:chislo_soch

nomer_soch=kk;

nc=setdiff(nodes,c(nomer_soch,:));% список не входящ в текущ сочет
nc(2,:)=0;

N=size(nc,2);% число вершин не входящ в текущ сочет

% прикрепляем вершины 

M3=M; 
M3(c(nomer_soch,:),:)=[];
%M3(:,nc(1,:))=[];
M3(:,nc(1,:))=100000;

[mina,nomer_prik]=min(M3,[],2);
nc(2,:)=nomer_prik';

if kk==1
   nc_vse=nc;
else
    nc_vse=vertcat(nc_vse,nc);
end

M2=M;
M2(nc(1,:),:)=[];
M2(:,nc(2,:))=[];

sumM2=sum(sum(M2,2));
R_c(nomer_soch,size(c,2)+1)=sumM2;

end

[min_R_c,nomer_stroki] = min(R_c(:,size(c,2)+1));

Pmedian=R_c(nomer_stroki,:);

SPmedian=cell(1,p);% Название пунктов для медиан и их прикрепленные
SPmedian(1,1:p)=raw(1,Pmedian(1:p)+1);
NPmedian(1,1:p)=Pmedian(1:p);
nc_min=nc;
nc_min(1:2,:)=nc_vse(2*nomer_stroki-1:2*nomer_stroki,:);

  for j=1:p 
     for i=1:size(nc_min,2)
         if (nc_min(2,i)==Pmedian(j))
             SPmedian(i+1,j)=raw(nc_min(1,i)+1,1);
             NPmedian(i+1,j)=nc_min(1,i);
         end
     end;
  end;

  [NN,MM]=size(SPmedian);
    
  for i=1:NN
      for j=1:MM
          if isempty(SPmedian(i,j))
            SPmedian(i,j)=' ';  
          end
      end
  end

% номера вершин p-медианы и минимальное суммарное от них расстояние
vse=zeros(size(nc_vse,1),size(nc_vse,2)+size(R_c,2));
vse(1:size(nc_vse,1),1:size(nc_vse,2))=nc_vse;
i=2:2:size(nc_vse,1);
vse(i,size(nc_vse,2)+1:size(nc_vse,2)+size(R_c,2))=R_c;
%четные строки номера фиксированных вершин (сочетания) и расстояния 
% нечетные прикрепленные вершины 
clear M M2 M3 MM N NN c chislo_soch i j...
    kk mina n nc nc_vse nodes nomer_prik R_c...
    nomer_soch nomer_stroki sumM2 vse;

[Koord_NP, NP] = shaperead('нас_пункты.shp');
[Koord_Dor, Dor] = shaperead('дороги.shp');

Ndor=length(Dor);
Nnp=length(NP);

for i=1:Nnp

    nazvNP(i).nazv=NP(i,1).TEXT_LABEL ;
     minj=Inf;    
 for j=1:Ndor

Ndorx=length(Koord_Dor(j,1).X);
xj=Koord_Dor(j,1).X(1,1:Ndorx-1);
yj=Koord_Dor(j,1).Y(1,1:Ndorx-1);
Rj=(xj-Koord_NP(i,1).X).^2+(yj-Koord_NP(i,1).Y).^2;
[zmin,pmin]=min(Rj);
if zmin<minj
   minj=zmin;
  nomer_dor(i,1)=i;
  nomer_dor(i,2)=j;
  
  
  nomer_dor(i,3)=pmin; % номер точки дороги
  nomer_dor(i,4)=minj;
  
end   
end
end
nazvanie_NP=char(nazvNP.nazv);
 s_nomer_dor=sortrows(nomer_dor,2);

 MR=zeros(Nnp,Nnp);
 for i=1:Nnp-1
     for j=i+1:Nnp
 
     if nomer_dor(i,2)==nomer_dor(j,2)
         MR(i,j)=nomer_dor(i,2);
     end        
     end
 end
         
 
Ndorx=length(Koord_Dor(6,1).X);
xj=Koord_Dor(6,1).X(1,1:Ndorx-1);
yj=Koord_Dor(6,1).Y(1,1:Ndorx-1);

%[xi,yi]=removeExtraNanSeparators(Koord_Dor(1,1).X,Koord_Dor(1,1).Y);
ko=241.1187;
 [tuva1, R] = geotiffread('tuva1.tif');
 mapshow(tuva1, R);
 geoshow('границы тувы1.shp');
 kk=166;
 
 ss1='{''COVER'',''грунтовая'',''Color'',[0 1 1]}';
 ss2='{''COVER'',''с покрытием'',''Color'',[1 0 1]}';
 %ss=strcat(ss1,',',ss2);

 symspecDor = makesymbolspec('Line',eval(ss1),eval(ss2));

 %  symspecDor = makesymbolspec('Line',...
% {'ROADS7_ID',166,'Color',[0 1 1]});     
  

geoshow('дороги.shp','SymbolSpec', symspecDor,...
 'DefaultColor',[0 0 1],...
  'DefaultLineWidth' ,2);
 
 symspecNP = makesymbolspec('Point', ...
   {'TEXT_LABEL','Кара-Холь', 'MarkerFaceColor', 'green','MarkerSize', 1,'Marker','o'}, ...
   {'TEXT_LABEL','Кызыл-Хая','MarkerFaceColor', 'magenta','MarkerSize', 1,'Marker','o'});



 
    set(gcf,'pos',[1 1 1920 1080]);
    

% lat = [Koord_NP.X]';
% lon = [Koord_NP.Y]';
% mat = char(NP.TEXT_LABEL);
% qrydata(gca,'City Data',{'NP','qrytest',lat,lon,mat});


[Nm,Mm]=size(SPmedian);

svet=jet(Mm);
%svet=hsv(Mm);


kk=2;
for k=1:Nnp
    for i=1:Nm
        for j=1:Mm
            if (strcmpi(NP(k).TEXT_LABEL,cell2mat(SPmedian(i,j)))==1)
            kk=kk+1; 
            XY(i,j,1)=Koord_NP(k).X;
            XY(i,j,2)=Koord_NP(k).Y;
            
            
   symspecNP.MarkerFaceColor(kk,1:3)={'TEXT_LABEL',NP(k).TEXT_LABEL,[svet(j,1) svet(j,2) svet(j,3)]};
   %symspecNP.MarkerFaceColor(kk,1:3)={'TEXT_LABEL',NP(k).TEXT_LABEL,[svetR(j) svetG(j) svetB(j)]};
      
         
         if i==1
         symspecNP.MarkerSize(kk,1:3)={'TEXT_LABEL',NP(k).TEXT_LABEL,14};
         symspecNP.Marker(kk,1:3)={'TEXT_LABEL',NP(k).TEXT_LABEL,'d'};
         
         else
         symspecNP.MarkerSize(kk,1:3)={'TEXT_LABEL',NP(k).TEXT_LABEL,8};
         symspecNP.Marker(kk,1:3)={'TEXT_LABEL',NP(k).TEXT_LABEL,'o'};
          end;    
             
        
            end
        end
    end
end


%xlswrite('Resultat1.xls',SPmedian);



    mapshow('нас_пункты.shp', 'SymbolSpec',symspecNP,...
   'DefaultMarkerEdgeColor', 'red',...
   'DefaultMarkerSize', 2, ...
   'DefaultMarker', 'o' ,...
   'DefaultMarkerFaceColor', 'blue');

hold on;
%axesm('MapProjection','utm','MapLatLimit',[46 54],'MapLonLimit',[70 98]);
%axesm('MapProjection','utm');

%,'MapLatLimit',[70 98],'MapLonLimit',[46 54])
%line([XY(1,1,1) ;XY(1,2,1)],[XY(1,1,2);XY(1,2,2) ],'r-');

% xk=find(XY(:,1,1)>0);
% yk=find(XY(:,1,2)>0);
% y=XY(xk,1,1);
% x=XY(yk,1,2);
% y1=y;
% y1(:,1)=y(1,1);
% x1=x;
% x1(:,1)=x(1,1);
%      plot(y,x,'--gs',...
%     'LineWidth',2,'MarkerEdgeColor','k',...
%                 'MarkerFaceColor','g',...
%                 'MarkerSize',2);
           


%plot(XY(1,1,1),XY(1,1,2),XY(1,2,1),XY(1,2,2),'r*');

% plot([46;48;50;52;54],[70;75;80;90;95],'g-');
% hold on;
[latt,lon] = track2(48,75,52,80);

% Plot both tracks.
%plot(98,52,'r*');

