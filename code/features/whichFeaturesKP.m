function fields = whichFeaturesKP(spectra)

fields = fieldnames(spectra);
do_keep = false(1,numel(fields));
for ii = 1:numel(fields)
    do_keep(ii) = size(spectra.(fields{ii}),2) == size(spectra.times,2) & ...
        ~strcmp(fields{ii},'times') & ~strcmp(fields{ii},'psd') & ~strcmp(fields{ii},'deriv');
end
fields(~do_keep) = [];
