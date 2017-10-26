function [files, dat,fileID] = compileStats

files = multiDirData('spikeSongStats*');
nSessions = numel(files);
fileID = cell(nSessions,1);
dat = cell(nSessions,1);
for ii = 1:nSessions
    [a,b]=strtok(files{ii},'\'); [a,b]=strtok(b,'\'); [a,b]=strtok(b,'-'); 
    a = strtok(b,'.'); 
    fileID{ii} = a(2:end);
    dat{ii} = load(files{ii});
end

nNeurons = zeros(nSessions,1);
baseRates = [];
colmap = jet(nSessions);
for ii = 1:nSessions
    nNeurons(ii) = numel(dat{ii}.stats);
    xx = [dat{ii}.stats.baseRate];
    yy = [dat{ii}.stats.juvieRate]; yySEM = [dat{ii}.stats.juvieRateSEM];
    zz = [dat{ii}.stats.adultRate]; zzSEM = [dat{ii}.stats.adultRateSEM];
    
    marker = '.';
    %if dat{ii}.stats
    figure(1)
    loglog(xx,zz,'.','MarkerEdgeColor', colmap(ii,:)); %  hold on;
    xlabel('Local baseline firing (Hz)');
    ylabel('Adult song firing (Hz)');
    hold on;

    figure(2)
    loglog(xx,yy,'.','MarkerEdgeColor', colmap(ii,:)); %  hold on;v
    xlabel('Local baseline firing (Hz)');
     ylabel('Juvenile song firing (Hz)');
    hold on;
end

IDleg = strrep(fileID,'_','-');
figure(1)
xl = xlim; yl = ylim; 
hold on; loglog(xl,xl,'k-'); hold off;
legend(IDleg);
figure(2)
xl = xlim; yl = ylim; 
hold on; loglog(xl,xl,'k-'); hold off;
legend(IDleg);

% unify the structure a little bit...
nTotNeurons = sum(nNeurons);
datCompiled = initEmptyStructArray(fieldnames(dat{1}.stats), nTotNeurons);
ctr = 1;
IDfield = cell(1,nTotNeurons);
ageField = cell(1,nTotNeurons);

recordMeta = readDataFile('RECORDINGDATA.txt');
keyfield = cell(1,numel(recordMeta));
for ii = 1:numel(recordMeta)
    keyfield{ii} = [recordMeta(ii).bird '_' recordMeta(ii).recording];
end

for ii = 1:nSessions
    assignRange = ctr:(ctr+nNeurons(ii)-1);
    datCompiled(assignRange) = dat{ii}.stats;
    [IDfield{assignRange}] = deal(fileID{ii});
    
    thisKey = find(strcmp(fileID{ii}, keyfield), 1)
    thisAge = recordMeta(thisKey).age;
    [ageField{assignRange}] = deal(thisAge);
    
    ctr = ctr + nNeurons(ii);
end
[datCompiled.ID] = IDfield{:};
[datCompiled.age] = ageField{:};

dat = datCompiled;


