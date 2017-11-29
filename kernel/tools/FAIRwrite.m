function [filename] = FAIRwrite(T,m,filename)
[p,f,e] = fileparts(filename);
filename = fullfile(p,[f,'.jpg']);
fprintf('write data to %s\n',filename)
imwrite(fliplr(uint8(reshape(T,m)))',filename,'jpeg');
%eval(['!open ',filename])