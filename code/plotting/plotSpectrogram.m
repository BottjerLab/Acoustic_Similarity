function hndl = plotSpectrogram(spec)
% formatting figure
spec.psd = spec.psd + eps;
surf(spec.times, spec.freqs, 10*log10(abs(spec.psd)),'EdgeColor','none');
view(0,90);

xlim([min(spec.times) max(spec.times)]); xlabel('Time (s)');
ylim([min(spec.freqs) max(spec.freqs)]); ylabel('Frequency (Hz)')

%colormap(hot);
caxis([-100,-25]);

%{
hold on;
      if isfield(spec,'centerFreq')
            plot(spec.times,spec.centerFreq, 'g-');
        end
        if isfield(spec,'harmonicPitch')
            plot(spec.times,spec.harmonicPitch, 'y-');
        end
  hold off;
  %}
%{
        hold on;
        if isfield(spec,'minCenterFreq')
            plot(spec.times,params.minCenterFreq,'b--');
        end
        hold off;
%}
hndl = gca;
end
