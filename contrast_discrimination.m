function behavior = contrast_discrimination(subjectID)
%%%%%%%%%%
% v1.0 20150710 Luke Hallum
%%%%%%%%%%

% Feedback tones dis/abled.
global MGL;

%%%
% Seed rand -- see 'help rand'.
%%%%%%%%%%
rand('twister',ceil(10000*rem(now,1)));
%%%%%%%%%%

%%%
% Design, display, target parameters. Anonymous functions.
%%%%%%%%%%
RECT_LEN_DEG = 5;
REF_CON = 0.5; % probably 0.5, meaning 75% gray on a 50% gray background
RECT_CON_MAX = 0.1; % max contrast difference between rectangles; this can't be too high
DURATION_PROBE_S = 2.0;
DURATION_BLANK1_S = 1.0;
WIDTH_STIM_DEG = 20;
WID_FIXATION_DEG = 0.2;
NUM_REPEATS = 5;
fnShuffle = @(x) x(randperm(length(x)));
vecDesignContrast = fnShuffle(RECT_CON_MAX*(repmat(-3.5:1:3.5,[1 NUM_REPEATS]) / 3.5));
vecDesignRotation = 2*(double(rand(size(vecDesignContrast)) > 0.5) - 0.5);
%%%%%%%%%%

%%%
% The format of matrix 'behavior'.
%%%%%%%%%%
IX_BEH_CON = 1;
IX_BEH_CORRECT = 2;
behavior = -1*ones(length(vecDesignContrast),2);
%%%%%%%%%%

%%%
% Get started...
%%%%%%%%%%
mglOpen();
mglVisualAngleCoordinates(57,[40 30]);
mglPolygon(WIDTH_STIM_DEG/2*[-1 1 1 -1],WIDTH_STIM_DEG/2*[-1 -1 1 1],0.5*[1 1 1]) % gray
mglPolygon(WID_FIXATION_DEG/2*[-1 1 1 -1],WID_FIXATION_DEG/2*[-1 -1 1 1],[1 0 0]) % fixation mark
mglFlush();
%%%%%%%%%%

%%%
% Experiment starts after a key press...
%%%%%%%%%%
KEYS_RESPONSE = [39 41]; % 'j' and 'k'
respKeys = 0*mglGetKeys();
disp('Press ''j'' or ''k'' to start experiment...')
while (sum(respKeys([KEYS_RESPONSE])) < 1), respKeys = mglGetKeys(); end
%%%%%%%%%%

for iiTrial = 1:length(vecDesignRotation)

  % Show blank, target, blank.
  mglPolygon(WIDTH_STIM_DEG/2*[-1 1 1 -1],WIDTH_STIM_DEG/2*[-1 -1 1 1],0.5*[1 1 1]) % gray
  mglPolygon(WID_FIXATION_DEG/2*[-1 1 1 -1],WID_FIXATION_DEG/2*[-1 -1 1 1],[1 0 0]) % fixation mark
  mglFlush();
  mglWaitSecs(DURATION_BLANK1_S);
  Rctangle(vecDesignRotation(iiTrial), RECT_LEN_DEG, REF_CON + vecDesignContrast(iiTrial)/2);
  Rctangle(-vecDesignRotation(iiTrial), RECT_LEN_DEG, REF_CON - vecDesignContrast(iiTrial)/2);
  mglFlush();
  mglWaitSecs(DURATION_PROBE_S);
  mglPolygon(WIDTH_STIM_DEG/2*[-1 1 1 -1],WIDTH_STIM_DEG/2*[-1 -1 1 1],0.5*[1 1 1]) % gray
  mglPolygon(WID_FIXATION_DEG/2*[-1 1 1 -1],WID_FIXATION_DEG/2*[-1 -1 1 1],[1 0 0]) % fixation mark
  mglFlush();

  % Get behavior.
  respKeys = 0*mglGetKeys();
  while (sum(respKeys([KEYS_RESPONSE])) < 1), respKeys = mglGetKeys(); end
  behavior(iiTrial,IX_BEH_CON) = vecDesignContrast(iiTrial);
  behavior(iiTrial,IX_BEH_CORRECT) = 0;
  if (respKeys(KEYS_RESPONSE(1)) == 1 & (vecDesignContrast(iiTrial) < 0)), behavior(iiTrial,IX_BEH_CORRECT) = 1; end
  if (respKeys(KEYS_RESPONSE(2)) == 1 & (vecDesignContrast(iiTrial) > 0)), behavior(iiTrial,IX_BEH_CORRECT) = 1; end

  % Give feedback.
  if(behavior(iiTrial,IX_BEH_CORRECT) > 0.5), mglPlaySound(find(strcmp(MGL.soundNames,'Blow'))); end

end

mglClose;

filename = sprintf('contrast_discrimination_%s_%s',subjectID,datestr(now,30));
save(filename,'behavior')

return;

%%%
% Some functions...
%%%%%%%%%%

function Rctangle(switchRotation, length_deg, contrast)
% Input 'length_deg' is that of the rectangle's longest side.

CANONICAL_VERTICAL_X = [0.2; -0.2; -0.2; 0.2];
CANONICAL_VERTICAL_Y = [1; 1; -1; -1];
ROTATION_DEG = 45;
% Inputs x and y are column vectors; output is n-by-2 (ie., horzcat of 2 column vectors).
fnRotation = @(x,y,rot_deg) transpose([cos(rot_deg/180*pi) -sin(rot_deg/180*pi); sin(rot_deg/180*pi) cos(rot_deg/180*pi)] * [x'; y']);

  this_rctangle = length_deg/2 * fnRotation(CANONICAL_VERTICAL_X, CANONICAL_VERTICAL_Y, switchRotation*ROTATION_DEG);
  mglPolygon(this_rctangle(:,1),this_rctangle(:,2),0.5 + contrast*[0.5 0.5 0.5]);

return;

