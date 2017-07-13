function behavior = illusory_contours(subjectID)
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
SCALE = 1.4; % scale of the canonical (ie., 1-deg diameter) Pacman/Varin figure
ECC_DEG = 5; % eccentricity from the stimulus center at which Pacmans/Varin figures appear
MOUTH_DEG = 90; % Pacmans/Varin figures have 'mouth' subtending this angle
DURATION_PROBE_S = 2.0;
DURATION_BLANK1_S = 1.0;
WIDTH_STIM_DEG = 20;
WID_FIXATION_DEG = 0.2;
NUM_REPEATS = 5;
vecDesignInducers = repmat(2:5, [1 NUM_REPEATS]);
vecDesignInducers = vecDesignInducers(randperm(length(vecDesignInducers)));
vecDesignRotation = 2*(double(rand(size(vecDesignInducers)) > 0.5) - 0.5);
%%%%%%%%%%

%%%
% The format of matrix 'behavior'.
%%%%%%%%%%
IX_BEH_NUM_IND = 1;
IX_BEH_ROT = 2;
IX_BEH_BUTTON = 3;
behavior = -1*ones(length(vecDesignRotation),3);
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
  Foursome(vecDesignRotation(iiTrial), SCALE, ECC_DEG, MOUTH_DEG, 0);
  Foursome(-vecDesignRotation(iiTrial), SCALE, ECC_DEG, MOUTH_DEG, vecDesignInducers(iiTrial));
  mglFlush();
  mglWaitSecs(DURATION_PROBE_S);
  mglPolygon(WIDTH_STIM_DEG/2*[-1 1 1 -1],WIDTH_STIM_DEG/2*[-1 -1 1 1],0.5*[1 1 1]) % gray
  mglPolygon(WID_FIXATION_DEG/2*[-1 1 1 -1],WID_FIXATION_DEG/2*[-1 -1 1 1],[1 0 0]) % fixation mark
  mglFlush();

  % Get behavior.
  respKeys = 0*mglGetKeys();
  while (sum(respKeys([KEYS_RESPONSE])) < 1), respKeys = mglGetKeys(); end
  behavior(iiTrial,IX_BEH_NUM_IND) = vecDesignInducers(iiTrial);
  behavior(iiTrial,IX_BEH_ROT) = vecDesignRotation(iiTrial);
  behavior(iiTrial,IX_BEH_BUTTON) = -1;
  if (respKeys(KEYS_RESPONSE(1)) == 1), behavior(iiTrial,IX_BEH_BUTTON) = 0; end
  if (respKeys(KEYS_RESPONSE(2)) == 1), behavior(iiTrial,IX_BEH_BUTTON) = 1; end

end

mglClose;

filename = sprintf('contrast_sensitivity_%s_%s',subjectID,datestr(now,30));
save(filename,'behavior')

return;

%%%
% Some functions...
%%%%%%%%%%

function Foursome(switchRotation, scale, eccentricityDeg, angleMouthDeg, numInducers)
% This function renders a Pacman/Varin figure foursome rotated either -45deg (switch = -1) or
% 45deg (switch = 1) from canonical ENWS.

CANONICAL_ENWS_X = [1; 0; -1; 0];
CANONICAL_ENWS_Y = [0; 1; 0; -1];
ROTATION_DEG = 22.5;
% Inputs x and y are column vectors; output is n-by-2 (ie., horzcat of 2 column vectors).
fnRotation = @(x,y,rot_deg) transpose([cos(rot_deg/180*pi) -sin(rot_deg/180*pi); sin(rot_deg/180*pi) cos(rot_deg/180*pi)] * [x'; y']);

  this_enws_xy = eccentricityDeg * fnRotation(CANONICAL_ENWS_X, CANONICAL_ENWS_Y, switchRotation*ROTATION_DEG);
  for (ii = 1:length(CANONICAL_ENWS_X))
    if (numInducers < 1), PacmanSingle(this_enws_xy(ii,1),this_enws_xy(ii,2),scale,angleMouthDeg,atan2(this_enws_xy(ii,2),this_enws_xy(ii,1))/pi*180); end
    if (numInducers > 0), VarinSingle(this_enws_xy(ii,1),this_enws_xy(ii,2),scale,angleMouthDeg,atan2(this_enws_xy(ii,2),this_enws_xy(ii,1))/pi*180,numInducers); end
  end

return;

function PacmanSingle(xPosDeg, yPosDeg, scale, angleMouthDeg, angleFacingDeg)
%
%

CONTRAST_PACMAN = 0.7;
NUM_SLICES = 36;
% Inputs x and y are column vectors; output is n-by-2 (ie., horzcat of 2 column vectors).
fnRotation = @(x,y,rot_deg) transpose([cos(rot_deg/180*pi) -sin(rot_deg/180*pi); sin(rot_deg/180*pi) cos(rot_deg/180*pi)] * [x'; y']);

  xy = scale*[cos(2*pi*transpose(linspace(0,1,NUM_SLICES+1))) sin(2*pi*transpose(linspace(0,1,NUM_SLICES+1)))];
  xy = xy - repmat([xPosDeg yPosDeg],[size(xy,1) 1]);
  mglPolygon(xy(:,1),xy(:,2),CONTRAST_PACMAN*[1 1 1]);
  xy = fnRotation(scale*[0; 1.1; 1.1], scale*[0; 1.1*tan(angleMouthDeg/2/180*pi); -1.1*tan(angleMouthDeg/2/180*pi)], angleFacingDeg);
  xy = xy - repmat([xPosDeg yPosDeg],[size(xy,1) 1]);
  mglPolygon(xy(:,1),xy(:,2),0.5*[1 1 1]);

return;

function VarinSingle(xPosDeg, yPosDeg, scale, angleMouthDeg, angleFacingDeg, numInducers)
%
%

NUM_SLICES = 36;
% Inputs x and y are column vectors; output is n-by-2 (ie., horzcat of 2 column vectors).
fnRotation = @(x,y,rot_deg) transpose([cos(rot_deg/180*pi) -sin(rot_deg/180*pi); sin(rot_deg/180*pi) cos(rot_deg/180*pi)] * [x'; y']);

  scales = scale * ([1:numInducers] / numInducers);
  for iiscale = 1:length(scales)
    xy = scales(iiscale)*[cos(2*pi*transpose(linspace(0,1,NUM_SLICES+1))) sin(2*pi*transpose(linspace(0,1,NUM_SLICES+1)))];
    xy = xy - repmat([xPosDeg yPosDeg],[size(xy,1) 1]);
    mglLines2(xy(1:(end-1),1),xy(1:(end-1),2),xy(2:end,1),xy(2:end,2),2,[1 1 1],1);
    xy = fnRotation(scales(iiscale)*[0; 1.1; 1.1], scales(iiscale)*[0; 1.1*tan(angleMouthDeg/2/180*pi); -1.1*tan(angleMouthDeg/2/180*pi)], angleFacingDeg);
    xy = xy - repmat([xPosDeg yPosDeg],[size(xy,1) 1]);
    mglPolygon(xy(:,1),xy(:,2),0.5*[1 1 1]);
  end

return;

