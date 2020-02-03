function volume4d = albira_readfiles(filename_hdr,path_hdr)
% read albira hdr and img raw files

cd(path_hdr)

info = albira_readhdr(filename_hdr);

[~,filename_img,~] = fileparts(filename_hdr);

volume4d = albira_readimg(filename_img, info);

end

function info = albira_readhdr(filename_hdr)
% read albira hdr text files
% extract basic information

fid = fopen(filename_hdr, 'r');
info = [];
tline = fgets(fid);
while ischar(tline)
    tline = fgets(fid);

    if strfind(tline, 'total_frames')
        C = strsplit(tline);
        info.total_frames = str2double(C{2});
    end
    if strfind(tline, 'x_dimension')
        C = strsplit(tline);
        info.x_dimension = str2double(C{2});
    end
    if strfind(tline, 'y_dimension')
        C = strsplit(tline);
        info.y_dimension = str2double(C{2});
    end
    if strfind(tline, 'z_dimension')
        C = strsplit(tline);
        info.z_dimension = str2double(C{2});
    end
    
end
fclose(fid);

end

function volume4d = albira_readimg(filename_img, info)
% read albira files
% use basic info

t = info.total_frames;
x = info.x_dimension;
y = info.y_dimension;
z = info.z_dimension;

fin = fopen(filename_img,'r');
I=fread(fin,x*y*z*t,'real*4','ieee-le');
volume4d = reshape(I,x,y,z,t);

[~,~,dim_z,dim_t] = size(volume4d);

for i=1:dim_z
    for j=1:dim_t        
        aux = rot90(rot90(rot90(volume4d(:,:,i,j))));
        volume4d(:,:,i,j) = fliplr(aux);
    end
end

fclose(fin);

end