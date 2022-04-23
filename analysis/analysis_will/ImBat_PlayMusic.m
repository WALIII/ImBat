function ImBat_PlayMusic(flightPaths);

% remember you need to pause the function in order to play in the
notecreate = @(frq,dur) sin(2*pi* [1:dur]/8192 * (440*2.^((frq-1)/12)));
notename = {'A' 'A#' 'B' 'C' 'C#' 'D' 'D#' 'E' 'F' 'F#' 'G' 'G#'};
song = {'A' 'A' 'E' 'E' 'F#' 'F#' 'E' 'E' 'D' 'D' 'C#' 'C#' 'B' 'B' 'A' 'A'};



scale2use = [1 4 6 7 8 11 13 16 18 19 20 23 25 23 20 19 18 16 13 11 8  7 6 4 1 4 6 7 8 11 13 16 18 19 20 23 25];
scale2use = cat(2, scale2use, 1:200);
%pentatonic = [1, 2, 3, 5, 6]
% sort flight paths by time
[iii idx2] = sort(flightPaths.flight_starts_idx);
flightPaths.id(idx2)

songidx = scale2use(flightPaths.id(idx2));
sinterval= diff(iii);
sinterval = [sinterval 1];


songnote = [];
for k1 = 1:length(songidx)
    songnote = [songnote; [notecreate(songidx(k1),1000)  zeros(1,75)]'];
    %songnote = [songnote; [notecreate(songidx(k1),1000)  zeros(1,round(sinterval(k1)/10))]'];

end
% 
player = audioplayer(songnote, 8192);

% break point here
play(player);

