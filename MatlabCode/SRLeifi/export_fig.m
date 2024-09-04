function [imageData, alpha] = export_fig(varargin) %#ok<*STRCL1,*DATST,*TNOW1>

% ======================
% Hjúpa export_fig með sama nafni í þessari möppu svo ég geti kallað á
% upprunalega fallið í undirmöppunni ExportFig með því að fara inn í þá
% möppu og svo til baka þegar það er búið.
%
% Þetta er gert því ég nenni ekki að fá öll þessi föll í aðalkeyrslumöppuna
% sem og að bæta henni sem matlab path
%
%
% A.t.h. ef export_fig klikkar þá þarf líklegast að keyra 
% cd ..
% handvirkt til að komast aftur í grunnmöppuna
% =======================


cd ./ExportFig/
[imageData, alpha] =  export_fig(varargin{:});
cd ..

