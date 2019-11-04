function [ Y,allnames] = readfromfolder(dataPATH,func,varargin)
%This function read all the data directly from the folder into variable Y
[iscell,variablename,npre] = process_options(varargin,'iscell','1','variablename',0,'npre',5);


allfiles = dir(dataPATH);


pres = 0;
for i = 1:npre
    s=allfiles(i).name;
    if regexp(s, '^\.')
        pres = pres+1;
    end
end

if ispc
    getname = @(i) [dataPATH ,'\' allfiles(i).name];
    getfilename = @(i) allfiles(i).name;
else
    getname = @(i) [dataPATH ,'/' allfiles(i).name];
    getfilename = @(i) allfiles(i).name;
end


n = length(allfiles);
allnames = cell(n-pres,1);
for i = (pres+1):n
    allnames{i-pres} = getfilename(i);
end

%Y1 = func(getname(pres+1));

Y = 0;
if ~isempty(func)
    if iscell
        Y = cell(n-pres,1);
        %    Y0 = func(getname(pres+1));
        for i = (pres+1):n
            Y1=func(getname(i));
            
            if ~variablename
                Y{i-pres} = Y1;
            else
                eval(['Y1 = Y1.',variablename,';'])
                Y{i-pres} = Y1;
            end
        end
        
    else
        Y0 = func(getname(pres+1));
        if ~variablename
            Y0 = Y0*0;
        else
            eval(['Y0 = Y0.',variablename,'*0;'])
        end
        Ysize = size(Y0);
        
        Y = repmat(Y0,[ones(1,length(size(Y0))-1),n-pres]);
        for i = (pres+1):n
            i
            ii = i-pres;
            Y1=func(getname(i));
            if ~variablename
                Y1 = Y1*0;
            else
                eval(['Y1 = Y1.',variablename,';'])
            end
            
            Y(:,:,((ii-1)*Ysize(end)+1):(ii*Ysize(end))) = Y1;
        end
        
        
    end
end

