function [t,Ca] = import_voistat(file)

filename = file;
delimiter = '\t';
startRow = 4;

formatSpec = '%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%s%*s%*s%*s%[^\n\r]';

fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);

% Allocate imported array to column variable names
t = dataArray{:, 4};
Ca = dataArray{:, 7};

end




