   
list1={'PSO', 'QPSO'};
% [indxPSO, tf]=listdlg('ListString',list1);

% d = dir;
% fn = {d.name};
[indxPSO, tf] = listdlg('PromptString',{'Bir algoritma se�iniz.',''},...
    'SelectionMode','single','ListString',list1);


list2={'Mevcut Durum Senaryosu', 'Do�al Gaz Yeni Yat�r�m Yasakl� Senaryo','Yeni K�m�r Santralleri Yasakl� Senaryo', ...
    'Yeni Fosil Santrali Yasakl� Senaryo', 'Karbon Vergisi Senaryosu(5$/tonCO2e)','Karbon Vergisi Senaryosu(10$/tonCO2e)',...
    'Karbon Vergisi Senaryosu(30$/tonCO2e)','Harici Maliyetli-MDS Senaryo', ...
    'Harici Maliyetli-Do�al Gaz Yeni Yat�r�m Yasakl� Senaryo','Harici Maliyetli-Yeni K�m�r Santralleri Yasakl� Senaryo', ...
    'Harici Maliyetli-Yeni Fosil Santrali Yasakl� Senaryo'};
% [indxScn, tf]=listdlg('ListString',list2);
[indxPSO, tf] = listdlg('PromptString',{'Bir senaryo se�iniz.',''},...
    'SelectionMode','single','ListString',list2);

list3={'Do�rusal Olan', 'Do�rusal Olmayan'};
%[indxLnr, tf]=listdlg('ListString',list3);
[indxPSO, tf] = listdlg('PromptString',{'Do�rusall�k �art�n� se�iniz.',''},...
    'SelectionMode','single','ListString',list3);
cc=5;
if indxPSO==1
    funcPSO_custom(indxPSO,indxScn,indxLnr);
elseif indxPSO==2    
    funcQPSO_custom(indxPSO,indxScn,indxLnr);
end

