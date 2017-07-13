function behavior = contrast_sensitivity(subjectID,signContrast)
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
% Design, display, target parameters.
%%%%%%%%%%
DURATION_PROBE_S = 0.2;
DURATION_BLANK1_S = 1.0;
WID_STIM_DEG = 20;
WID_TARGET_DEG = 2;
WID_FIXATION_DEG = 0.2;
ECC_TARGET_DEG = 4;
NUM_REPEATS = 5;
MAX_CONTRAST = 0.15;
MIN_CONTRAST = 0.01;
NUM_CONTRAST = 10;
vecDesignCon = signContrast*repmat(linspace(MIN_CONTRAST,MAX_CONTRAST,NUM_CONTRAST),[1 NUM_REPEATS]);
vecDesignCon = vecDesignCon(randperm(length(vecDesignCon)));
vecDesignLR = 2*(double(rand(1,length(vecDesignCon)) > 0.5) - 0.5);
%%%%%%%%%%

%%%
% The format of matrix 'behavior'.
%%%%%%%%%%
IX_BEH_CON = 1;
IX_BEH_LR = 2;
IX_BEH_BUTTON = 3;
IX_BEH_CORRECT = 4;
behavior = -1*ones(length(vecDesignCon),4);
%%%%%%%%%%

%%%
% Get started...
%%%%%%%%%%
mglOpen();
mglVisualAngleCoordinates(57,[40 30]);
mglPolygon(WID_STIM_DEG/2*[-1 1 1 -1],WID_STIM_DEG/2*[-1 -1 1 1],0.5*[1 1 1]) % gray
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

for iiTrial = 1:length(vecDesignCon)

  % Show blank, target, blank.
  mglPolygon(WID_STIM_DEG/2*[-1 1 1 -1],WID_STIM_DEG/2*[-1 -1 1 1],0.5*[1 1 1]) % gray
  mglPolygon(WID_FIXATION_DEG/2*[-1 1 1 -1],WID_FIXATION_DEG/2*[-1 -1 1 1],[1 0 0]) % fixation mark
  mglFlush();
  mglWaitSecs(DURATION_BLANK1_S);
  X_JIT = 2*(rand-0.5) + vecDesignLR(iiTrial)*ECC_TARGET_DEG; %% Target position is jittered in x and y.
  Y_JIT = 2*(rand-0.5);                                       %% Include target position in this jitter.
  mglPolygon(WID_STIM_DEG/2*[-1 1 1 -1],WID_STIM_DEG/2*[-1 -1 1 1],0.5*[1 1 1]) % gray
  mglPolygon(WID_FIXATION_DEG/2*[-1 1 1 -1],WID_FIXATION_DEG/2*[-1 -1 1 1],[1 0 0]) % fixation mark
  mglPolygon(WID_TARGET_DEG/2*[-1 1 1 -1]+X_JIT,WID_TARGET_DEG/2*[-1 -1 1 1]+Y_JIT,vecDesignCon(iiTrial)/2+[0.5 0.5 0.5]) % target
  mglFlush();
  mglWaitSecs(DURATION_PROBE_S);
  mglPolygon(WID_STIM_DEG/2*[-1 1 1 -1],WID_STIM_DEG/2*[-1 -1 1 1],0.5*[1 1 1]) % gray
  mglPolygon(WID_FIXATION_DEG/2*[-1 1 1 -1],WID_FIXATION_DEG/2*[-1 -1 1 1],[1 0 0]) % fixation mark
  mglFlush();

  % Get behavior.
  respKeys = 0*mglGetKeys();
  while (sum(respKeys([KEYS_RESPONSE])) < 1), respKeys = mglGetKeys(); end
  behavior(iiTrial,IX_BEH_CON) = vecDesignCon(iiTrial);
  behavior(iiTrial,IX_BEH_LR) = vecDesignLR(iiTrial);
  behavior(iiTrial,IX_BEH_BUTTON) = -1;
  if (respKeys(KEYS_RESPONSE(1)) == 1), behavior(iiTrial,IX_BEH_BUTTON) = 0; end
  if (respKeys(KEYS_RESPONSE(2)) == 1), behavior(iiTrial,IX_BEH_BUTTON) = 1; end
  behavior(iiTrial,IX_BEH_CORRECT) = 0;
  if ((vecDesignLR(iiTrial) == -1) & (behavior(iiTrial,IX_BEH_BUTTON) == 0)), behavior(iiTrial,IX_BEH_CORRECT) = 1; end
  if ((vecDesignLR(iiTrial) == 1) & (behavior(iiTrial,IX_BEH_BUTTON) == 1)), behavior(iiTrial,IX_BEH_CORRECT) = 1; end

  % Give feedback.
  if(behavior(iiTrial,IX_BEH_CORRECT) > 0.5), mglPlaySound(find(strcmp(MGL.soundNames,'Blow'))); end

end

mglClose;

filename = sprintf('contrast_sensitivity_%s_%s',subjectID,datestr(now,30));
save(filename,'behavior')

return;

