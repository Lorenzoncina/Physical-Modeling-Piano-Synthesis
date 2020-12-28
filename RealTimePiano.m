clear;

%visualizza o meno il plot della stringa
plot = true;

%impostazione GUI con html
fig = uifigure;
fig.Position = [10 10 800 400]; 
h = uihtml(fig);
h.Position = [0 0 800 400];
h.HTMLSource = fullfile(pwd,'gui2.html');
%h.DataChangedFcn = @(src,event)disp(h.Data);
%midiNote = h.Data;

%gestione virtual keyboard
%h.DataChangedFcn = @(src,event)PianoSoundSynthesis(midi(h.Data+1));
h.DataChangedFcn = @(src,event)copyData(jsondecode(h.Data()));

%function needed to pass parameters to PianoSoundSynthesis function
function [nota,riverbero] = copyData(txt)
    %conversione midi/frequency
    a = 440; 
    midi = zeros(128);
    %notare che la nota midi 0 si trova alla posizione 1, la nota midi 1 alla
    %posizione 2 ecc..
    for x = 1:128
      midi(x) = (a/32)*(2^(((x-1)-9)/12));
    end
    nota = str2num(txt.midi);
    riverbero = txt.reverbType;
    disp(nota);
    disp(riverbero);
    PianoSoundSynthesis(midi(nota+1), riverbero);
end


%gestione midi keybord collegata
% 
% device = mididevice('CASIO USB-MIDI');
% while true
%     %gestione tastiera midi
%     receivedMessages = midireceive(device);
%     for i=1:numel(receivedMessages)
%         msg = receivedMessages(i);
%         if isNoteOn(msg)
%             notaMidi = msg.Note;
%             fundamentalFrequency = midi(notaMidi+1);
%             PianoSoundSynthesis(fundamentalFrequency, plot);
%         end
%     end
% end
%  
% function yes = isNoteOn(msg)
%     yes = msg.Type == midimsgtype.NoteOn ...
%         && msg.Velocity > 0;
% end



    



