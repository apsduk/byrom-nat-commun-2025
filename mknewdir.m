function fname = mknewdir(fname)
% This function creates a new folder with the name fname.
% (c) Copyright Dan Byrom and Alexander Darlington 2024.

% --- Determine file seperator --------------------------------------------
if ispc
    sep = '\';
elseif isunix
    sep = '/';
elseif ismac
    sep = '/';
else
    error('Unsupported OS');
end

% ----- Make sure fname is string -----------------------------------------
fname = char(fname);

% ----- List current files named fname ------------------------------------
list = dir([pwd,sep,fname,'*']);

% ----- Identify if isdir -------------------------------------------------
isfname = zeros(length(list),1);
isdir   = zeros(length(list),1);
for i = 1:length(list)
    
    % ----- Get filename --------------------------------------------------
    txt = strsplit(list(i).name,'_runid_');
    if length(txt) == 1
        txt = '';
    else
        txt = txt(end-1);
    end
    
    % ----- Save logical tests --------------------------------------------
    isfname(i) = strcmp(txt, fname);
    isdir(i) = list(i).isdir;
    
end

% ----- Number of existing dir --------------------------------------------
if isempty(list)
    N = 0;
else
    N = sum(isdir(logical(isfname)));
end

% ----- Build new dir name ------------------------------------------------
fname = [fname,'_runid_',num2str(N+1)];

% ----- Make new directory ------------------------------------------------
mkdir(fname);

end
