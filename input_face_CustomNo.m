   
list1={'PSO', 'QPSO'};
% [indxPSO, tf]=listdlg('ListString',list1);

% d = dir;
% fn = {d.name};
[indxPSO, tf] = listdlg('PromptString',{'Bir algoritma seçiniz.',''},...
    'SelectionMode','single','ListString',list1);


list2={'Mevcut Durum Senaryosu', 'Doðal Gaz Yeni Yatýrým Yasaklý Senaryo','Yeni Kömür Santralleri Yasaklý Senaryo', ...
    'Yeni Fosil Santrali Yasaklý Senaryo', 'Karbon Vergisi Senaryosu(5$/tonCO2e)','Karbon Vergisi Senaryosu(10$/tonCO2e)',...
    'Karbon Vergisi Senaryosu(30$/tonCO2e)','Harici Maliyetli-MDS Senaryo', ...
    'Harici Maliyetli-Doðal Gaz Yeni Yatýrým Yasaklý Senaryo','Harici Maliyetli-Yeni Kömür Santralleri Yasaklý Senaryo', ...
    'Harici Maliyetli-Yeni Fosil Santrali Yasaklý Senaryo'};
% [indxScn, tf]=listdlg('ListString',list2);
[indxPSO, tf] = listdlg('PromptString',{'Bir senaryo seçiniz.',''},...
    'SelectionMode','single','ListString',list2);

list3={'Doðrusal Olan', 'Doðrusal Olmayan'};
%[indxLnr, tf]=listdlg('ListString',list3);
[indxPSO, tf] = listdlg('PromptString',{'Doðrusallýk þartýný seçiniz.',''},...
    'SelectionMode','single','ListString',list3);
cc=5;
if indxPSO==1
    funcPSO_custom(indxPSO,indxScn,indxLnr);
elseif indxPSO==2    
    funcQPSO_custom(indxPSO,indxScn,indxLnr);
end

