function data = loadAsciiData(file)
%  data = loadAsciiData(file)
%     Load data from the AMP ascii writer.
% 
%  This function loads the results from AMP's ascii writer into matlab.
%  file is the filename to load.  The results are returned as an array
%  of structures.  The format for each entry is:
%     data.type = {'Vector','Matrix'}
%     data.name = "name of the vector or matrix"
%     data.data = the actual data.  For matrices this will be in sparse form

fid = fopen(file);
data = struct([]);
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    i = length(data)+1;
    if strfind(tline,'Vector:')==1
        % Loading a vector
        data(i).type = 'Vector';
        index = find(tline=='"');
        data(i).name = tline(index(1)+1:index(2)-1);
        N = str2double(tline(index(2)+1:length(tline)));
        data(i).data = zeros(N,1);
        for j = 1:N
            tline = fgetl(fid);
            data(i).data(j) = str2double(tline);
        end
    elseif strfind(tline,'Matrix:')==1
        % Load a matrix
        data(i).type = 'Matrix';
        index = find(tline=='"');
        data(i).name = tline(index(1)+1:index(2)-1);
        N = str2num(tline(index(2)+1:length(tline))); %#ok<ST2NM>
        data(i).data = sparse(N(1),N(2));
        tline = fgetl(fid);
        while ~isempty(tline)
            x = str2num(tline); %#ok<ST2NM>
            data(i).data(x(1)+1,x(2)+1) = x(3);
            tline = fgetl(fid);
        end
    elseif isempty(tline)
        % Empty line, ignore
    else
        error('Unknown data detected');
    end
end
fclose(fid);


