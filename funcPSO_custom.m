function funcPSO_custom(indxPSO,indxScn,indxLnr)
tic
bireyselHizKatsayisi =2;
suruHizKatsayisi     =2;
genelHizKatsayisi    =0.8;
bireySayisi          =500;
time                 =16;
iterasyonSayisi      =1000;
plantType            =7;
cFact                =[0.25 0.15 0.85 0.75 0.85 0.85 0.85];
payda                =[100 10 50 50 300 300 300];
optQ                 =[8 1 2 3 10 10 10];
%                     Wind PV Bio Geo Lig NGcc NGbsc

eff                  =[1 1 0.39 0.15 0.38 0.56 0.39];
invCostPrm           =[173.106 242.3488 311.5913 255.9315 216.12 69.832 60.3095];
fCostPrm             =[21.38 30.081	100	100	37.818 30.568 34];
vCostPrm             =[0 0 0 0 3 6.5816 6.5816];
vCost_hPrm           =[0.2487 2.85 4.976 4.976 10.4628 4.976 7.8683];
emissionPrm          =[0 0 111.84 138.89 95.52 53.06 53.06];
talepUstSinir        =[52.77 57.23 62.05 67.31 72.99 79.18 85.87 93.13 101 109.54 118.81 128.86 139.75 151.56 164.37 178.26];
K_Vergisi_5          =[0 0 0.005 0.005 0.005  0.005 0.005];
K_Vergisi_10         =[0 0 0.01  0.01  0.01   0.01  0.01];
K_Vergisi_30         =[0 0 0.03 0.03 0.03  0.03 0.03]; 
suruDegerItr=1000000000*ones(1,iterasyonSayisi);
altSinir=0;
B=100000*ones(1,16);
clc;
objIt=1000000*ones(iterasyonSayisi);
    if (indxScn==1)||(indxScn==8)
        plantCap             =[50.765 4.7304 10.722 9.4608 1000 1000 1000];
    elseif (indxScn==2)||(indxScn==9)
        plantCap             =[50.765 4.7304 10.722 9.4608 1000 7.7117 0];
    elseif (indxScn==3)||(indxScn==10)
        plantCap             =[50.765 4.7304 10.722 9.4608 53.8945 1000 1000];
    elseif (indxScn==4)||(indxScn==11)
        plantCap             =[50.765 4.7304 10.722 9.4608 53.8945 7.7117 0];
    elseif (indxScn>4)&&(indxScn<8)
        plantCap             =[50.765 4.7304 10.722 9.4608 1000 1000 1000];
    end
ctp=0;
tprm=zeros(1,time);
tprm(1,1:time)=(1/1.05).^(0:time-1);
for i=1:time
%     tprm(1,i)=(1/1.05)^(i-1);
     if sum(plantCap(1,:))<talepUstSinir(1,i)
         ctp=ctp+1;
     end
 end
if ctp==0
    fprintf('Plan Olurludur');
else
    fprintf('Plan Olursuzdur');
end

if ctp==0
for oRun=1:1
%list={'BAU Senaryo', 'Doðal Gaz Yeni Yatýrým Yasaklý Senaryo','Yeni Kömür Santralleri Yasaklý Senaryo', ...
   % 'Yeni Fosil Santrali Yasaklý Senaryo', 'Karbon Vergisi Senaryosu','Harici Maliyetli-BAU Senaryo', ...
   % 'Harici Maliyetli-Doðal Gaz Yeni Yatýrým Yasaklý Senaryo','Harici Maliyetli-Yeni Kömür Santralleri Yasaklý Senaryo', ...
   % 'Harici Maliyetli-Yeni Fosil Santrali Yasaklý Senaryo'};
suru =zeros(plantType, time, bireySayisi);
obj  =zeros(bireySayisi,1);
for k=1:bireySayisi
    if (indxScn==1)||(indxScn==8)
        plantCap             =[50.765 4.7304 10.722 9.4608 1000 1000 1000];
    elseif (indxScn==2)||(indxScn==9)
        plantCap             =[50.765 4.7304 10.722 9.4608 1000 7.7117 0];
    elseif (indxScn==3)||(indxScn==10)
        plantCap             =[50.765 4.7304 10.722 9.4608 53.8945 1000 1000];
    elseif (indxScn==4)||(indxScn==11)
        plantCap             =[50.765 4.7304 10.722 9.4608 53.8945 7.7117 0];
    elseif (indxScn>4)&&(indxScn<8)
        plantCap             =[50.765 4.7304 10.722 9.4608 1000 1000 1000];
    end
    
    
    
    
    
%     suru(1:plantType,1:time,k) =rand*(min(plantCap(1,1:plantType), talepUstSinir(1,1:time))-altSinir)+altSinir;
    for t=1:time
        for j=1:plantType
%             suru(j,t,k) =unifrnd(altSinir,min(plantCap(1,j), talepUstSinir(1,t)));
            suru(j,t,k) =rand*(min(plantCap(1,j), talepUstSinir(1,t))-altSinir)+altSinir;
        end
    end
    plantAvgCap=zeros(plantType,time);
    plantAvgCap(1,1)=6.60222;
    plantAvgCap(2,1)=0.2167;
    plantAvgCap(3,1)=0;  
    plantAvgCap(4,1)=2.790936;         
    plantAvgCap(5,1)=53.8945;
    plantAvgCap(6,1)=7.71117;
    plantAvgCap(7,1)=0; 
    
    for j=1:time
        totalutilized=sum(suru(:,j,k));
        while talepUstSinir(1,j)-totalutilized>0
              R=fix(rand*(plantType)+1);              
              %R = randperm(plantType,1); 
              % remove those four numbers from A
              if talepUstSinir(1,j)-totalutilized>plantCap(1,R)
                 suru(R,j,k)=suru(R,j,k)+plantCap(1,R);  
                 totalutilized=totalutilized+plantCap(1,R);
              else
                 suru(R,j,k)=suru(R,j,k)+talepUstSinir(1,j)-totalutilized;                 
                 totalutilized=talepUstSinir(1,j);
              end   
        end
        boyut=zeros(1,plantType);
        while totalutilized>talepUstSinir(1,j) 
            ctr=0;
%             clear('boyut');            
            for m=1:plantType
                if suru(m,j,k)>0
                   ctr=ctr+1;
                   boyut(1,ctr)=m;    
                end
                if m==plantType
                    boyut=boyut(1:ctr);
                end
            end
              randIdcs = randperm(length(boyut),1);            
              R = boyut(1,randIdcs);
              % remove those four numbers from A
              if totalutilized-talepUstSinir(1,j)>suru(R,j,k)
                 totalutilized=totalutilized-suru(R,j,k);
                 suru(R,j,k)=0; 
              else
                 suru(R,j,k)=suru(R,j,k)-(totalutilized-talepUstSinir(1,j));                 
                 totalutilized=talepUstSinir(1,j);
              end
%               clear('boyut');
        end 
        
        
        
    end        
                    
    %Her santralin maks kap 1 olacaktý 1'i aþan var mý diye bakýyoruz  

    t=1;
    while t<=time 
        i=1;
        while i<=plantType 
        if suru(i,t,k)>plantCap(1,i) 
             fark(1,t)=suru(i,t,k)-plantCap(1,i);
             suru(i,t,k)=plantCap(1,i);                
                boyut=zeros(1,plantType);
                while fark(1,t)>0
                    ctr=0;
%                     clear('boyut');
                    for m=1:plantType
                        if (plantCap(1,m)-suru(m,t,k))>0.000001
                            ctr=ctr+1;
                            boyut(1,ctr)=m;                    
                        end
                        if m==plantType
                            boyut=boyut(1:ctr);
                        end
                    end
                    %randIdcs = randperm(length(boyut),1);     
                    randIdcs=fix(rand*(length(boyut))+1);
                    R = boyut(1,randIdcs);
                    if fark(1,t)<=(plantCap(1,R)-suru(R,t,k))
                        suru(R,t,k)=suru(R,t,k)+fark(1,t);
                        fark(1,t)=0;
                    else
                        fark(1,t)=fark(1,t)-(plantCap(1,R)-suru(R,t,k));
                        suru(R,t,k)=plantCap(1,R);
                    end                    
                    
                end
     
                i=0;
                end           
                i=i+1;
        end
            for mm=1:plantType
                if t>1
                    plantAvgCap(mm,t)=plantAvgCap(mm,t-1);
                end
                if plantAvgCap(mm,t)<suru(mm,t,k)
                    plantAvgCap(mm,t)=suru(mm,t,k);             
                end
            end
    t=t+1;
end   
       i=k;          
            obj(i)=0  ;    

            if (indxScn<5)&&(indxLnr==1)
                for n2=1:time
                    for n=1:plantType
                        obj(i)=obj(i)+(plantAvgCap(n,n2)/31.536/cFact(1,n))*(fCostPrm(1,n)+invCostPrm(1,n))*tprm(1,n2);                          
                        obj(i)=obj(i)+(vCostPrm(1,n)/eff(1,n))*suru(n,n2,k)*tprm(1,n2); 
                   end   
                end
            elseif (indxScn<5)&&(indxLnr==2)
                for n2=1:time
                    for n=1:plantType
                        obj(i)=obj(i)+(plantAvgCap(n,n2)/31.536/cFact(1,n))*(((1+(((optQ(1,n)-(plantAvgCap(n,n2)/31.536/cFact(1,n)))^2)/payda(1,n)))*fCostPrm(1,n))+invCostPrm(1,n))*tprm(1,n2);                          
                        obj(i)=obj(i)+(vCostPrm(1,n)/eff(1,n))*suru(n,n2,k)*tprm(1,n2); 
                   end   
                end
            elseif (indxScn==5)&&(indxLnr==1)
                for n2=1:time
                    for n=1:plantType
                        obj(i)=obj(i)+(plantAvgCap(n,n2)/31.536/cFact(1,n))*(fCostPrm(1,n)+invCostPrm(1,n))*tprm(1,n2);                          
                        obj(i)=obj(i)+(vCostPrm(1,n)/eff(1,n)+emissionPrm(1,n)*K_Vergisi_5(1,n)/eff(1,n))*suru(n,n2,k)*tprm(1,n2); 
                   end   
                end
            elseif (indxScn==6)&&(indxLnr==1)
                for n2=1:time
                    for n=1:plantType
                        obj(i)=obj(i)+(plantAvgCap(n,n2)/31.536/cFact(1,n))*(fCostPrm(1,n)+invCostPrm(1,n))*tprm(1,n2);                          
                        obj(i)=obj(i)+(vCostPrm(1,n)/eff(1,n)+emissionPrm(1,n)*K_Vergisi_10(1,n)/eff(1,n))*suru(n,n2,k)*tprm(1,n2); 
                   end   
                end
            elseif (indxScn==7)&&(indxLnr==1)
                for n2=1:time
                    for n=1:plantType
                        obj(i)=obj(i)+(plantAvgCap(n,n2)/31.536/cFact(1,n))*(fCostPrm(1,n)+invCostPrm(1,n))*tprm(1,n2);                          
                        obj(i)=obj(i)+(vCostPrm(1,n)/eff(1,n)+emissionPrm(1,n)*K_Vergisi_30(1,n)/eff(1,n))*suru(n,n2,k)*tprm(1,n2); 
                   end   
                end
            elseif (indxScn==5)&&(indxLnr==2)
                for n2=1:time
                    for n=1:plantType
                        obj(i)=obj(i)+(plantAvgCap(n,n2)/31.536/cFact(1,n))*(((1+(((optQ(1,n)-(plantAvgCap(n,n2)/31.536/cFact(1,n)))^2)/payda(1,n)))*fCostPrm(1,n))+invCostPrm(1,n))*tprm(1,n2);                          
                        obj(i)=obj(i)+(vCostPrm(1,n)/eff(1,n)+emissionPrm(1,n)*K_Vergisi_5(1,n)/eff(1,n))*suru(n,n2,k)*tprm(1,n2);
                   end   
                end   
            elseif (indxScn==6)&&(indxLnr==2)
                for n2=1:time
                    for n=1:plantType
                        obj(i)=obj(i)+(plantAvgCap(n,n2)/31.536/cFact(1,n))*(((1+(((optQ(1,n)-(plantAvgCap(n,n2)/31.536/cFact(1,n)))^2)/payda(1,n)))*fCostPrm(1,n))+invCostPrm(1,n))*tprm(1,n2);                          
                        obj(i)=obj(i)+(vCostPrm(1,n)/eff(1,n)+emissionPrm(1,n)*K_Vergisi_10(1,n)/eff(1,n))*suru(n,n2,k)*tprm(1,n2);
                   end   
                end   
            elseif (indxScn==7)&&(indxLnr==2)
                for n2=1:time
                    for n=1:plantType
                        obj(i)=obj(i)+(plantAvgCap(n,n2)/31.536/cFact(1,n))*(((1+(((optQ(1,n)-(plantAvgCap(n,n2)/31.536/cFact(1,n)))^2)/payda(1,n)))*fCostPrm(1,n))+invCostPrm(1,n))*tprm(1,n2);                          
                        obj(i)=obj(i)+(vCostPrm(1,n)/eff(1,n)+emissionPrm(1,n)*K_Vergisi_30(1,n)/eff(1,n))*suru(n,n2,k)*tprm(1,n2);
                   end   
                end 
            elseif (indxScn>7)&&(indxLnr==1)
                for n2=1:time
                    for n=1:plantType
                        obj(i)=obj(i)+(plantAvgCap(n,n2)/31.536/cFact(1,n))*(fCostPrm(1,n)+invCostPrm(1,n))*tprm(1,n2);                          
                        obj(i)=obj(i)+(vCostPrm(1,n)/eff(1,n)+vCost_hPrm(1,n))*suru(n,n2,k)*tprm(1,n2);
                   end   
                end
            elseif (indxScn>7)&&(indxLnr==2)
                for n2=1:time
                    for n=1:plantType
                        obj(i)=obj(i)+(plantAvgCap(n,n2)/31.536/cFact(1,n))*(((1+(((optQ(1,n)-(plantAvgCap(n,n2)/31.536/cFact(1,n)))^2)/payda(1,n)))*fCostPrm(1,n))+invCostPrm(1,n))*tprm(1,n2);                          
                        obj(i)=obj(i)+(vCostPrm(1,n)/eff(1,n)+vCost_hPrm(1,n))*suru(n,n2,k)*tprm(1,n2);
                   end   
                end
            end
            
end

surudekiEnIyiBireyinYeri=zeros(plantType,time);
suruEnIyiDegeri=min(obj);
idx=find(suruEnIyiDegeri==obj);
surudekiEnIyiBireyinYeri=suru(:,:,idx);
hiz=zeros(plantType,t,bireySayisi);
bireyselEnIyiPozisyon=suru;
bireyselEnIyiDeger=min(obj);
objIt(1)=suruEnIyiDegeri;


iterasyon=1;
while(iterasyon <= iterasyonSayisi)
    hiz=zeros(plantType,time,bireySayisi); 
    for k=1:bireySayisi  
        hiz(:,:,k) = genelHizKatsayisi * hiz(:,:,k) + bireyselHizKatsayisi * rand * (bireyselEnIyiPozisyon(:,:,k)- suru(:,:,k));
        hiz(:,:,k) = hiz(:,:,k)                     + suruHizKatsayisi     * rand * (surudekiEnIyiBireyinYeri(:,:) - suru(:,:,k));    
    if (indxScn==1)||(indxScn==8)
        plantCap             =[50.765 4.7304 10.722 9.4608 1000 1000 1000];
    elseif (indxScn==2)||(indxScn==9)
        plantCap             =[50.765 4.7304 10.722 9.4608 1000 7.7117 0];
    elseif (indxScn==3)||(indxScn==10)
        plantCap             =[50.765 4.7304 10.722 9.4608 53.8945 1000 1000];
    elseif (indxScn==4)||(indxScn==11)
        plantCap             =[50.765 4.7304 10.722 9.4608 53.8945 7.7117 0];
    elseif (indxScn>4)&&(indxScn<8)
        plantCap             =[50.765 4.7304 10.722 9.4608 1000 1000 1000];
    end
    
    
%         hiz(:,:,k)=(-1*suru(:,:,k))*(hiz(:,:,k)<(-1*suru(:,:,k)))+hiz(:,:,k)*(hiz(:,:,k)>=(-1*suru(:,:,k)));        
%         hiz(:,:,k)=(min(plantCap(1,i),talepUstSinir(1,j))-suru(i,j,k))*(hiz(i,j,k)>min(plantCap(1,i),talepUstSinir(1,j))-suru(i,j,k))
        for j=1:time
            for i=1:plantType
                if hiz(i,j,k)<-1*suru(i,j,k)
                   hiz(i,j,k)=-1*suru(i,j,k);
                end
%                 hiz(i,j,k)=(-1*suru(i,j,k))*(hiz(i,j,k)<(-1*suru(i,j,k)))+hiz(i,j,k)*(hiz(i,j,k)>=(-1*suru(i,j,k)));
%                 hiz(:,:,k)=(min(plantCap(1,i),talepUstSinir(1,j))-suru(i,j,k))*(hiz(i,j,k)>min(plantCap(1,i),talepUstSinir(1,j))-suru(i,j,k))+hiz(i,j,k)*((hiz(i,j,k)<=min(plantCap(1,i),talepUstSinir(1,j))-suru(i,j,k))&&(hiz(i,j,k)>0));
                if hiz(i,j,k)>min(plantCap(1,i),talepUstSinir(1,j))-suru(i,j,k)
                   hiz(i,j,k)=min(plantCap(1,i),talepUstSinir(1,j))-suru(i,j,k);
                end
                suru(i,j,k) =suru(i,j,k)+hiz(i,j,k);
            end
        end 

     for j=1:time
            totalutilized=sum(suru(:,j,k));
            if totalutilized>talepUstSinir(1,j)
                boyut=zeros(1,plantType);
              while totalutilized-talepUstSinir(1,j)>0
                ctr=0;
%                     clear('boyut');
                    for n=1:plantType
                        if suru(n,j,k)>0
                            ctr=ctr+1;
                            boyut(1,ctr)=n;
                        end
                        if n==plantType
                            boyut=boyut(1:ctr);
                        end
                    end                
                    %R=boyut(randperm(length(boyut),1));
                    R= boyut(1,fix(rand*(length(boyut))+1));
                    if suru(R,j,k)>=totalutilized-talepUstSinir(1,j)
                        suru(R,j,k)=suru(R,j,k)-(totalutilized-talepUstSinir(1,j));
                        totalutilized=talepUstSinir(1,j);
                    else
                        totalutilized=totalutilized-suru(R,j,k);
                        suru(R,j,k)=0;
                    end
%                     clear('boyut');
                end   
            end
            if talepUstSinir(1,j)-totalutilized>0
                while talepUstSinir(1,j)-totalutilized>0
                    %R=randperm(plantType,1);
                    R= fix(rand*(plantType)+1);
                      if talepUstSinir(1,j)-totalutilized>plantCap(1,R)-suru(R,j,k);
                         totalutilized=totalutilized+plantCap(1,R)-suru(R,j,k);
                         suru(R,j,k)=plantCap(1,R);  
                      else
                         suru(R,j,k)=suru(R,j,k)+talepUstSinir(1,j)-totalutilized;                 
                         totalutilized=talepUstSinir(1,j); 
                      end 
                end
             end
     end   

   
mutasyonPrb0 =unifrnd(0,1);
if mutasyonPrb0<=0.1  

for t=1:time 
        mutasyonPrb =unifrnd(0,1);
        if mutasyonPrb<=0.1
            r = randi([1 plantType],1,1);
            m_pi=unifrnd(0,1);
            sinpi=sin(pi*m_pi);
            suru(r,t,k)=sinpi*min(talepUstSinir(1,t),plantCap(1,r));
%             suru(r,t,k)=unifrnd(0,min(talepUstSinir(1,t),plantCap(1,r)));
            if suru(r,t,k)<0
                suru(r,t,k)=0;
            end
        end
        totalutilized=sum(suru(:,t,k));
        while talepUstSinir(1,t)-totalutilized>0
              R = randperm(plantType,1);               
              totalutilized=totalutilized+plantCap(1,R)-suru(R,t,k);
              suru(R,t,k)=plantCap(1,R);  
              if talepUstSinir(1,t)-totalutilized<=plantCap(1,R)-suru(R,t,k)
                 suru(R,t,k)=suru(R,t,k)+talepUstSinir(1,t)-totalutilized;                 
                 totalutilized=talepUstSinir(1,t);
              end   
        end
        while totalutilized>talepUstSinir(1,t) 
            ctr=0;
            clear('boyut');
            for m=1:plantType
                if suru(m,t,k)>0
                   ctr=ctr+1;
                   boyut(1,ctr)=m;                    
                end
            end
              randIdcs = randperm(length(boyut),1);            
              R = boyut(1,randIdcs);
              % remove those four numbers from A
              totalutilized=totalutilized-suru(R,t,k);
              suru(R,t,k)=0; 
              if totalutilized-talepUstSinir(1,t)<=suru(R,t,k)
                 suru(R,t,k)=suru(R,t,k)-(totalutilized-talepUstSinir(1,t));                 
                 totalutilized=talepUstSinir(1,t);
              end
              clear('boyut');
        end 
        
end
end
        
        plantAvgCap(1,1)=6.60222;
        plantAvgCap(2,1)=0.2167;
        plantAvgCap(3,1)=0;  
        plantAvgCap(4,1)=2.790936;         
        plantAvgCap(5,1)=53.8945;
        plantAvgCap(6,1)=7.71117;
        plantAvgCap(7,1)=0;     
        t=1;
        while t<=time             
            for mm=1:plantType 
                if t>1
                    plantAvgCap(mm,t)=plantAvgCap(mm,t-1);
                end
                if plantAvgCap(mm,t)<suru(mm,t,k)
                    plantAvgCap(mm,t)=suru(mm,t,k);             
                end
            end
        t=t+1;
        end
        
       i=k;           
       obj(i)=0  ;          
            if (indxScn<5)&&(indxLnr==1)
                for n2=1:time
                    for n=1:plantType                    
                        obj(i)=obj(i)+(plantAvgCap(n,n2)/31.536/cFact(1,n))*(fCostPrm(1,n)+invCostPrm(1,n))*tprm(1,n2);                          
                        obj(i)=obj(i)+(vCostPrm(1,n)/eff(1,n))*suru(n,n2,k)*tprm(1,n2); 
                   end   
                end
            elseif (indxScn<5)&&(indxLnr==2)
                for n2=1:time
                    for n=1:plantType
                        obj(i)=obj(i)+(plantAvgCap(n,n2)/31.536/cFact(1,n))*(((1+(((optQ(1,n)-(plantAvgCap(n,n2)/31.536/cFact(1,n)))^2)/payda(1,n)))*fCostPrm(1,n))+invCostPrm(1,n))*tprm(1,n2);                          
                        obj(i)=obj(i)+(vCostPrm(1,n)/eff(1,n))*suru(n,n2,k)*tprm(1,n2); 
                   end   
                end
            elseif (indxScn==5)&&(indxLnr==1)
                for n2=1:time
                    for n=1:plantType
                        obj(i)=obj(i)+(plantAvgCap(n,n2)/31.536/cFact(1,n))*(fCostPrm(1,n)+invCostPrm(1,n))*tprm(1,n2);                          
                        obj(i)=obj(i)+(vCostPrm(1,n)/eff(1,n)+emissionPrm(1,n)*K_Vergisi_5(1,n)/eff(1,n))*suru(n,n2,k)*tprm(1,n2); 
                   end   
                end
            elseif (indxScn==6)&&(indxLnr==1)
                for n2=1:time
                    for n=1:plantType
                        obj(i)=obj(i)+(plantAvgCap(n,n2)/31.536/cFact(1,n))*(fCostPrm(1,n)+invCostPrm(1,n))*tprm(1,n2);                          
                        obj(i)=obj(i)+(vCostPrm(1,n)/eff(1,n)+emissionPrm(1,n)*K_Vergisi_10(1,n)/eff(1,n))*suru(n,n2,k)*tprm(1,n2); 
                   end   
                end
            elseif (indxScn==7)&&(indxLnr==1)
                for n2=1:time
                    for n=1:plantType
                        obj(i)=obj(i)+(plantAvgCap(n,n2)/31.536/cFact(1,n))*(fCostPrm(1,n)+invCostPrm(1,n))*tprm(1,n2);                          
                        obj(i)=obj(i)+(vCostPrm(1,n)/eff(1,n)+emissionPrm(1,n)*K_Vergisi_30(1,n)/eff(1,n))*suru(n,n2,k)*tprm(1,n2); 
                   end   
                end
            elseif (indxScn==5)&&(indxLnr==2)
                for n2=1:time
                    for n=1:plantType
                        obj(i)=obj(i)+(plantAvgCap(n,n2)/31.536/cFact(1,n))*(((1+(((optQ(1,n)-(plantAvgCap(n,n2)/31.536/cFact(1,n)))^2)/payda(1,n)))*fCostPrm(1,n))+invCostPrm(1,n))*tprm(1,n2);                          
                        obj(i)=obj(i)+(vCostPrm(1,n)/eff(1,n)+emissionPrm(1,n)*K_Vergisi_5(1,n)/eff(1,n))*suru(n,n2,k)*tprm(1,n2);
                   end   
                end   
            elseif (indxScn==6)&&(indxLnr==2)
                for n2=1:time
                    for n=1:plantType
                        obj(i)=obj(i)+(plantAvgCap(n,n2)/31.536/cFact(1,n))*(((1+(((optQ(1,n)-(plantAvgCap(n,n2)/31.536/cFact(1,n)))^2)/payda(1,n)))*fCostPrm(1,n))+invCostPrm(1,n))*tprm(1,n2);                          
                        obj(i)=obj(i)+(vCostPrm(1,n)/eff(1,n)+emissionPrm(1,n)*K_Vergisi_10(1,n)/eff(1,n))*suru(n,n2,k)*tprm(1,n2);
                   end   
                end   
            elseif (indxScn==7)&&(indxLnr==2)
                for n2=1:time
                    for n=1:plantType
                        obj(i)=obj(i)+(plantAvgCap(n,n2)/31.536/cFact(1,n))*(((1+(((optQ(1,n)-(plantAvgCap(n,n2)/31.536/cFact(1,n)))^2)/payda(1,n)))*fCostPrm(1,n))+invCostPrm(1,n))*tprm(1,n2);                          
                        obj(i)=obj(i)+(vCostPrm(1,n)/eff(1,n)+emissionPrm(1,n)*K_Vergisi_30(1,n)/eff(1,n))*suru(n,n2,k)*tprm(1,n2);
                   end   
                end 
            elseif (indxScn>7)&&(indxLnr==1)
                for n2=1:time
                    for n=1:plantType
                        obj(i)=obj(i)+(plantAvgCap(n,n2)/31.536/cFact(1,n))*(fCostPrm(1,n)+invCostPrm(1,n))*tprm(1,n2);                          
                        obj(i)=obj(i)+(vCostPrm(1,n)/eff(1,n)+vCost_hPrm(1,n))*suru(n,n2,k)*tprm(1,n2);
                   end   
                end
            elseif (indxScn>7)&&(indxLnr==2)
                for n2=1:time
                    for n=1:plantType
                        obj(i)=obj(i)+(plantAvgCap(n,n2)/31.536/cFact(1,n))*(((1+(((optQ(1,n)-(plantAvgCap(n,n2)/31.536/cFact(1,n)))^2)/payda(1,n)))*fCostPrm(1,n))+invCostPrm(1,n))*tprm(1,n2);                          
                        obj(i)=obj(i)+(vCostPrm(1,n)/eff(1,n)+vCost_hPrm(1,n))*suru(n,n2,k)*tprm(1,n2);
                   end   
                end
            end
               
            if (obj(k) < suruEnIyiDegeri)
                suruEnIyiDegeri = obj(k);
                surudekiEnIyiBireyinYeri= suru(:,:,k);
            end   
%             annualEms=zeros(1,time);
            
            if (obj(k) < bireyselEnIyiDeger)
                bireyselEnIyiDeger = obj(k);
                bireyselEnIyiPozisyon=suru;
            end
%      fprintf('k=%f,\n',k);        
    end
    
    
%     suruDegerItr(1,iterasyon)=suruEnIyiDegeri;
%     if iterasyon>100
%         if suruDegerItr(1,iterasyon)>0.99*suruDegerItr(1,iterasyon-100)
%             iterasyonSayisi=iterasyon;
%             objIt=objIt(1:iterasyon);
%         end
%     end
    if iterasyon==iterasyonSayisi
        
       [ annualEms ] = emissionDetector(surudekiEnIyiBireyinYeri,eff,emissionPrm,plantType,time);
       tInterval=toc;
%        sheet = oRun;
%     B = cell2mat( myCell )
       A=suru(:,:,1);
       for k1=2:bireySayisi
    %         B{1,1} =strcat('Birey No:', num2str(k1));
            F=ones(1,16)*tInterval;
            A=vertcat(A,B,suru(:,:,k1));
            D=round(zeros(1,16)+suruEnIyiDegeri,8);
            E=ones(1,16)*iterasyon;
       end
           %format long;
           C=[ annualEms ];           
           A=vertcat(A,B,surudekiEnIyiBireyinYeri,B,C,B,D,B,E,B,F);           
           senaryoAdi=strcat(num2str(indxPSO),'_Senaryo_');
           senaryoAdi=strcat(senaryoAdi,num2str(indxScn));
           senaryoAdi=strcat(senaryoAdi,'_');
           senaryoAdi=strcat(senaryoAdi,num2str(indxLnr));
           senaryoAdi=strcat(senaryoAdi,'_');
           senaryoAdi=strcat(senaryoAdi,num2str(oRun));
           senaryoAdi=strcat(senaryoAdi,'.csv');
           dlmwrite(senaryoAdi,A,'precision','%.4f');
    end
    
    iterasyon=iterasyon+1;
    objIt(iterasyon)=suruEnIyiDegeri;
    fprintf('run no=%f, iterasyon no= %f\n',oRun,iterasyon);

genelHizKatsayisi=genelHizKatsayisi*0.99;    

%     
%         for j=1:time
%             if sum(suru(:,j,k))-talepUstSinir(1,j)>0.00001
%                 cb=0;
%             end
%         end
end % while
end
end
%toc
%     if ctp==0
%         plot(objIt); hold on;
%     end
end