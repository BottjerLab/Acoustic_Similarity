function clipAP = playSound(clip, sampleRate, doBlock)
% function CLIPAP = PLAYSOUND(CLIP, SAMPLERATE, DOBLOCK)
% clip is the vector of sound samples
% sampleRate is the sample rate of the sound file
% if doBlock is true, then this function waits for the sound to play before
% returning
% if doBlock is false, then this f
% the soundsc function also works, but we use this to control blocking
% behavior
% 
%   WARNING: This function will bluescreen on some machines -
%   make sure your sound card works with MATLAB
if nargin == 1 && strcmp('stop',clip)
end
if nargin < 3, doBlock = false; end;
if isnan(sampleRate), 
    error('playSound:badArgs', 'Please give a non-NaN sample rate');
end
try
    clipAP = audioplayer(clip, sampleRate);
    if ~doBlock
        played = false;
        while ~played
            try
                play(clipAP);
                played = true;
            catch ME
                if strcmp(ME.identifier,...
                        'MATLAB:audiovideo:audioplayer:DeviceError')
                    %fprintf('hitting warning...\n');
                    %keyboard
                    pause(0.25)
                    %sound(0,8192) % try clearing the sound buffer?
                end
            end
        end
    else
        playblocking(clipAP);
    end
catch ME
    clear clipAP;
    warning('Sound aborted, cleaning up...');
end
end
