
%-----------------------------------------
% set up ImSAnE path
%-----------------------------------------

% entries of matlab search path
pathentries = regexp(path, pathsep, 'split');

% find ones containing ImSAnE and remove from path
disp('removing old ImSAnE entries from path');
for i = 1:length(pathentries)
    if ~isempty(regexp(pathentries{i}, 'ImSAnE','once'))
        rmpath(pathentries{i});
    end
end

% path of current script is new ImSAnE directory
[ImSAnEpath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(ImSAnEpath));
disp('added ImSAnE directory containing setup to path');

% storing the path in the startup directory
upath = userpath;
upath = upath(1:end-1);
savepath(fullfile(upath, 'pathdef.m'));

%  add path to settings for easy retrieval
setpref('ImSAnE', 'path', ImSAnEpath);

%%
%-----------------------------------------
% settings
%-----------------------------------------

% detail of the level of output messages 
% admittedly not very well implemented
% 1: function names
% 2: function details
% 3: map call
msgLevel = 2;
setpref('ImSAnE', 'msgLevel', msgLevel);

% store figures to check quality of surface fits
fitQualityPlots =  1;
setpref('ImSAnE', 'fitQualityPlots', fitQualityPlots);

