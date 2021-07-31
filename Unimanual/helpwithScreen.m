Screen()
Usage:

% Activate compatibility mode: Try to behave like the old MacOS-9 Psychtoolbox:
oldEnableFlag=Screen('Preference', 'EmulateOldPTB', [enableFlag]);

% Open or close a window or texture:
[windowPtr,rect]=Screen('OpenWindow',windowPtrOrScreenNumber [,color] [,rect] [,pixelSize] [,numberOfBuffers] [,stereomode] [,multisample][,imagingmode][,specialFlags][,clientRect][,fbOverrideRect][,vrrParams=[]]);
[windowPtr,rect]=Screen('OpenOffscreenWindow',windowPtrOrScreenNumber [,color] [,rect] [,pixelSize] [,specialFlags] [,multiSample]);
textureIndex=Screen('MakeTexture', WindowIndex, imageMatrix [, optimizeForDrawAngle=0] [, specialFlags=0] [, floatprecision] [, textureOrientation=0] [, textureShader=0]);
oldParams = Screen('PanelFitter', windowPtr [, newParams]);
Screen('Close', [windowOrTextureIndex or list of textureIndices/offscreenWindowIndices]);
Screen('CloseAll');

%  Draw lines and solids like QuickDraw and DirectX (OS 9 and Windows):
currentbuffer = Screen('SelectStereoDrawBuffer', windowPtr [, bufferid] [, param1]);
Screen('DrawLine', windowPtr [,color], fromH, fromV, toH, toV [,penWidth]);
Screen('DrawArc',windowPtr,[color],[rect],startAngle,arcAngle)
Screen('FrameArc',windowPtr,[color],[rect],startAngle,arcAngle[,penWidth] [,penHeight] [,penMode])
Screen('FillArc',windowPtr,[color],[rect],startAngle,arcAngle)
Screen('FillRect', windowPtr [,color] [,rect] );
Screen('FrameRect', windowPtr [,color] [,rect] [,penWidth]);
Screen('FillOval', windowPtr [,color] [,rect] [,perfectUpToMaxDiameter]);
Screen('FrameOval', windowPtr [,color] [,rect] [,penWidth] [,penHeight] [,penMode]);
Screen('FramePoly', windowPtr [,color], pointList [,penWidth]);
Screen('FillPoly', windowPtr [,color], pointList [, isConvex]);

% New OpenGL functions for OS X:
Screen('glPoint', windowPtr, color, x, y [,size]);
Screen('gluDisk', windowPtr, color, x, y [,size]);
[minSmoothPointSize, maxSmoothPointSize, minAliasedPointSize, maxAliasedPointSize] = Screen('DrawDots', windowPtr, xy [,size] [,color] [,center] [,dot_type] [,lenient]);
[minSmoothLineWidth, maxSmoothLineWidth, minAliasedLineWidth, maxAliasedLineWidth] = Screen('DrawLines', windowPtr, xy [,width] [,colors] [,center] [,smooth] [,lenient]);
[sourceFactorOld, destinationFactorOld, colorMaskOld]=Screen('BlendFunction', windowIndex, [sourceFactorNew], [destinationFactorNew], [colorMaskNew]);

% Draw Text in windows
textModes = Screen('TextModes');
oldCopyMode=Screen('TextMode', windowPtr [,textMode]);
oldTextSize=Screen('TextSize', windowPtr [,textSize]);
oldStyle=Screen('TextStyle', windowPtr [,style]);
[oldFontName,oldFontNumber,oldTextStyle]=Screen('TextFont', windowPtr [,fontNameOrNumber][,textStyle]);
[normBoundsRect, offsetBoundsRect, textHeight, xAdvance] = Screen('TextBounds', windowPtr, text [,x] [,y] [,yPositionIsBaseline] [,swapTextDirection]);
[newX, newY, textHeight]=Screen('DrawText', windowPtr, text [,x] [,y] [,color] [,backgroundColor] [,yPositionIsBaseline] [,swapTextDirection]);
oldTextColor=Screen('TextColor', windowPtr [,colorVector]);
oldTextBackgroundColor=Screen('TextBackgroundColor', windowPtr [,colorVector]);
oldMatrix = Screen('TextTransform', windowPtr [, newMatrix]);

% Copy an image, very quickly, between textures, offscreen windows and onscreen windows.
[resident [texidresident]] = Screen('PreloadTextures', windowPtr [, texids]);
Screen('DrawTexture', windowPointer, texturePointer [,sourceRect] [,destinationRect] [,rotationAngle] [, filterMode] [, globalAlpha] [, modulateColor] [, textureShader] [, specialFlags] [, auxParameters]);
Screen('DrawTextures', windowPointer, texturePointer(s) [, sourceRect(s)] [, destinationRect(s)] [, rotationAngle(s)] [, filterMode(s)] [, globalAlpha(s)] [, modulateColor(s)] [, textureShader] [, specialFlags] [, auxParameters]);
Screen('CopyWindow', srcWindowPtr, dstWindowPtr, [srcRect], [dstRect], [copyMode])

% Copy an image, slowly, between matrices and windows :
imageArray=Screen('GetImage', windowPtr [,rect] [,bufferName] [,floatprecision=0] [,nrchannels=3])
Screen('PutImage', windowPtr, imageArray [,rect]);

% Synchronize with the window's screen (on-screen only):
[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed Beampos] = Screen('Flip', windowPtr [, when] [, dontclear] [, dontsync] [, multiflip]);
[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed Beampos] = Screen('AsyncFlipBegin', windowPtr [, when] [, dontclear] [, dontsync] [, multiflip]);
[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed Beampos] = Screen('AsyncFlipEnd', windowPtr);
[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed Beampos] = Screen('AsyncFlipCheckEnd', windowPtr);
[VBLTimestamp StimulusOnsetTime swapCertainTime] = Screen('WaitUntilAsyncFlipCertain', windowPtr);
[info] = Screen('GetFlipInfo', windowPtr [, infoType=0] [, auxArg1]);
[telapsed] = Screen('DrawingFinished', windowPtr [, dontclear] [, sync]);
framesSinceLastWait = Screen('WaitBlanking', windowPtr [, waitFrames]);

% Load color lookup table of the window's screen (on-screen only):
[gammatable, dacbits, reallutsize] = Screen('ReadNormalizedGammaTable', windowPtrOrScreenNumber [, physicalDisplay]);
[oldtable, success] = Screen('LoadNormalizedGammaTable', windowPtrOrScreenNumber, table [, loadOnNextFlip][, physicalDisplay][, ignoreErrors]);
oldclut = Screen('LoadCLUT', windowPtrOrScreenNumber [, clut] [, startEntry=0] [, bits=8]);

% Get (and set) information about a window or screen:
screenNumbers=Screen('Screens' [, physicalDisplays]);
windowPtrs=Screen('Windows');
kind=Screen(windowPtr, 'WindowKind');
isOffscreen=Screen(windowPtr,'IsOffscreen');
hz=Screen('FrameRate', windowPtrOrScreenNumber [, mode] [, reqFrameRate]);
hz=Screen('NominalFrameRate', windowPtrOrScreenNumber [, mode] [, reqFrameRate]);
[ monitorFlipInterval nrValidSamples stddev ]=Screen('GetFlipInterval', windowPtr [, nrSamples] [, stddev] [, timeout]);
screenNumber=Screen('WindowScreenNumber', windowPtr);
rect=Screen('Rect', windowPtrOrScreenNumber [, realFBSize=0]);
pixelSize=Screen('PixelSize', windowPtrOrScreenNumber);
pixelSizes=Screen('PixelSizes', windowPtrOrScreenNumber);
[width, height]=Screen('WindowSize', windowPointerOrScreenNumber [, realFBSize=0]);
[width, height]=Screen('DisplaySize', ScreenNumber);
[oldmaximumvalue, oldclampcolors, oldapplyToDoubleInputMakeTexture] = Screen('ColorRange', windowPtr [, maximumvalue][, clampcolors][, applyToDoubleInputMakeTexture]);
info = Screen('GetWindowInfo', windowPtr [, infoType=0] [, auxArg1]);
resolutions=Screen('Resolutions', screenNumber);
oldResolution=Screen('Resolution', screenNumber [, newwidth] [, newheight] [, newHz] [, newPixelSize] [, specialMode]);
oldSettings = Screen('ConfigureDisplay', setting, screenNumber, outputId [, newwidth][, newheight][, newHz][, newX][, newY]);
Screen('ConstrainCursor', windowIndex, addConstraint [, rect]);

% Get/set details of environment, computer, and video card (i.e. screen):
struct=Screen('Version');
comp=Screen('Computer');
oldBool=Screen('Preference', 'IgnoreCase' [,bool]);
tick0Secs=Screen('Preference', 'Tick0Secs', tick0Secs);
psychTableVersion=Screen('Preference', 'PsychTableVersion');
mexFunctionName=Screen('Preference', 'PsychTableCreator');
proc=Screen('Preference', 'Process');
Screen('Preference','SkipSyncTests', skipTest);
Screen('Preference','VisualDebugLevel', level (valid values between 0 and 5));
Screen('Preference', 'ConserveVRAM', mode (valid values between 0 and 3));
Screen('Preference', 'Enable3DGraphics', [enableFlag]);

% Helper functions.  Don't call these directly, use eponymous wrappers:
[x, y, buttonVector, hasKbFocus, valuators]= Screen('GetMouseHelper', numButtons [, screenNumber][, mouseIndex]);
Screen('HideCursorHelper', windowPntr [, mouseIndex]);
Screen('ShowCursorHelper', windowPntr [, cursorshapeid][, mouseIndex]);
Screen('SetMouseHelper', windowPntrOrScreenNumber, x, y [, mouseIndex][, detachFromMouse]);
Screen('SetMouseHelper', windowPntrOrScreenNumber, x, y [, mouseIndex][, detachFromMouse]);

% Internal testing of Screen
timeList= Screen('GetTimelist');
Screen('ClearTimelist');
Screen('Preference','DebugMakeTexture', enableDebugging);

% Movie and multimedia playback functions:
[ moviePtr [duration] [fps] [width] [height] [count] [aspectRatio] [hdrStaticMetaData]]=Screen('OpenMovie', windowPtr, moviefile [, async=0] [, preloadSecs=1] [, specialFlags1=0][, pixelFormat=4][, maxNumberThreads=-1][, movieOptions]);
Screen('CloseMovie' [, moviePtr=all]);
[ texturePtr [timeindex]]=Screen('GetMovieImage', windowPtr, moviePtr, [waitForImage], [fortimeindex], [specialFlags = 0] [, specialFlags2 = 0]);
[droppedframes] = Screen('PlayMovie', moviePtr, rate, [loop], [soundvolume]);
timeindex = Screen('GetMovieTimeIndex', moviePtr);
[oldtimeindex] = Screen('SetMovieTimeIndex', moviePtr, timeindex [, indexIsFrames=0]);
moviePtr = Screen('CreateMovie', windowPtr, movieFile [, width][, height][, frameRate=30][, movieOptions][, numChannels=4][, bitdepth=8]);
Screen('FinalizeMovie', moviePtr);
Screen('AddFrameToMovie', windowPtr [,rect] [,bufferName] [,moviePtr=0] [,frameduration=1]);
Screen('AddAudioBufferToMovie', moviePtr, audioBuffer);
[imageArray, format, errorMsg, auxInfo] = Screen('ReadHDRImage', filename [, errorMode=0]);

% Video capture functions:
devices = Screen('VideoCaptureDevices' [, engineId]);
videoPtr =Screen('OpenVideoCapture', windowPtr [, deviceIndex][, roirectangle][, pixeldepth][, numbuffers][, allowfallback][, targetmoviename][, recordingflags][, captureEngineType][, bitdepth=8]);
Screen('CloseVideoCapture', capturePtr);
[fps starttime] = Screen('StartVideoCapture', capturePtr [, captureRateFPS] [, dropframes=0] [, startAt]);
droppedframes = Screen('StopVideoCapture', capturePtr [, discardFrames=1]);
[ texturePtr [capturetimestamp] [droppedcount] [average_intensityOrRawImageMatrix]]=Screen('GetCapturedImage', windowPtr, capturePtr [, waitForImage=1] [,oldTexture] [,specialmode] [,targetmemptr]);
oldvalue = Screen('SetVideoCaptureParameter', capturePtr, 'parameterName' [, value]);

% Low level direct access to OpenGL-API functions:
% Online info for each function available by opening a terminal window
% and typing 'man Functionname' + Enter.

Screen('glPushMatrix', windowPtr);
Screen('glPopMatrix', windowPtr);
Screen('glLoadIdentity', windowPtr);
Screen('glTranslate', windowPtr, tx, ty [, tz]);
Screen('glScale', windowPtr, sx, sy [, sz]);
Screen('glRotate', windowPtr, angle, [rx = 0], [ry = 0] ,[rz = 1]);

% Support for 3D graphics rendering and for interfacing with external OpenGL code:
Screen('Preference', 'Enable3DGraphics', [enableFlag]);  % Enable 3D gfx support.
Screen('BeginOpenGL', windowPtr [, sharecontext]);  % Prepare window for external OpenGL drawing.
Screen('EndOpenGL', windowPtr);  % Finish external OpenGL drawing.
[targetwindow, IsOpenGLRendering] = Screen('GetOpenGLDrawMode');
[textureHandle rect] = Screen('SetOpenGLTextureFromMemPointer', windowPtr, textureHandle, imagePtr, width, height, depth [, upsidedown][, target][, glinternalformat][, gltype][, extdataformat][, specialFlags]);
[textureHandle rect] = Screen('SetOpenGLTexture', windowPtr, textureHandle, glTexid, target [, glWidth][, glHeight][, glDepth][, textureShader][, specialFlags]);
[ gltexid gltextarget texcoord_u texcoord_v ] =Screen('GetOpenGLTexture', windowPtr, textureHandle [, x][, y]);

% Support for plugins and for builtin high performance image processing pipeline:
[ret1, ret2, ...] = Screen('HookFunction', windowPtr, 'Subcommand', 'HookName', arg1, arg2, ...);
proxyPtr = Screen('OpenProxy', windowPtr [, imagingmode]);
transtexid = Screen('TransformTexture', sourceTexture, transformProxyPtr [, sourceTexture2][, targetTexture][, specialFlags]);