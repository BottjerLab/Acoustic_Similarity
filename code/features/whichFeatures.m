function fields = whichfeatures(spectra)

fields = fieldnames(spectra);
do_keep = false(1,numel(fields));
for ii = 1:numel(fields)
    do_keep(ii) = all(size(spectra.(fields{ii})) == size(spectra.times)) & ...
        ~strcmp(fields{ii},'times');
end
fields(~do_keep) = [];
