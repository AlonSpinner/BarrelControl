function Sfcn_DrawBeam(block)
setup(block);
end
function setup(block) %runs at t=0 i/o definitions
block.SetSimViewingDevice(true);

%dialog parameters
block.NumDialogPrms = 1;

%register number of ports
block.NumInputPorts = 3;
block.NumOutputPorts = 0;

%setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;

%Register the properties of the input ports

%Enable
block.InputPort(1).Complexity     ='Real';
block.InputPort(1).DataTypeId     =-1;
block.InputPort(1).Dimensions     =1;
block.InputPort(1).SamplingMode   ='Sample';

%q
block.InputPort(2).Complexity     ='Real';
block.InputPort(2).DataTypeId     =-1;
block.InputPort(2).Dimensions     = 8;
block.InputPort(2).SamplingMode   ='Sample';

%Time
block.InputPort(3).Complexity     ='Real';
block.InputPort(3).DataTypeId     =-1;
block.InputPort(3).Dimensions     =1;
block.InputPort(3).SamplingMode   ='Sample';

dtDraw=1/100;
block.SampleTimes = [dtDraw 0]; %[discrete time, offset]

%specify block simStateCompliace
block.SimStateCompliance = 'HasNoSimState';

%register functions
block.RegBlockMethod('InitializeConditions',    @InitializeConditions);
block.RegBlockMethod('Start',                   @Start);
block.RegBlockMethod('Terminate',               @Terminate);
block.RegBlockMethod('Outputs',                 @Outputs);
block.RegBlockMethod('CheckParameters',         @CheckPrms);
block.RegBlockMethod('ProcessParameters',       @ProcessPrms);
end
function ProcessPrms(block) %runs on every dt (Wasnt checked!)
  block.AutoUpdateRuntimePrms;
end
function InitializeConditions(block) %runs on t=0 and when susbystem is enabled
Enable=block.InputPort(1).Data(1);
if ~Enable, return, end

%check if figute exists and valid. if not - reset it
UserData=get(gcbh,'UserData');
if isempty(UserData) %first time simulation is activated
     SetupFigAndUserData(block);
elseif ~ishghandle(UserData.Figure) %figure was deleted
    SetupFigAndUserData(block);
else %figure exists, just clear it and start a new
    SetupFigAndUserData(block,UserData.Figure); %reset figure
end
end
function Outputs(block) %runs on every dt
UserData=get(gcbh,'UserData');
if ~ishghandle(UserData.Figure)
     UserData=SetupFigAndUserData(block); %set figure to a new start
end

%------Draw Measured Data
%Beam
q=block.InputPort(2).Data;
L=block.DialogPrm(1).Data;
[w_mm,xiVec_mm]=Compute_w(q,L,100);
UserData.hBeam.XData=xiVec_mm;
UserData.hBeam.YData=w_mm;

UserData.hforce.YData=[w_mm(1),20];

%Spring
[x,y]=UserData.oSpr.getSpr([0,w_mm(1)],[0,-30]);
UserData.hSpr.XData=x;
UserData.hSpr.YData=y;

%Update time text
Time=block.InputPort(3).Data(1);
UserData.hTime.String=sprintf('Time %g[s]',Time);

drawnow limitrate
end
%% Auxiliary functions
function UserData=SetupFigAndUserData(block,varargin)
q=block.InputPort(2).Data;
L=block.DialogPrm(1).Data;

if nargin<2 %figure was not provided in input
    %Create figure
    FigureName='OnlyPhysics';
    Fig = figure(...
        'Name',              FigureName,...
        'NumberTitle',        'off',...
        'IntegerHandle',     'off',...
        'Color',             [1,1,1],...
        'MenuBar',           'figure',...
        'ToolBar',           'auto',...
        'HandleVisibility',   'callback',...
        'Resize',            'on',...
        'visible',           'on');
    
    %Create Axes
    Ax=axes(Fig);
    hold(Ax,'on'); grid(Ax,'on');
    xlim(Ax,1e3*[-1,L]); %in mm
    ylim(Ax,50*[-1,1]); %in mm
%     axis(Ax,'equal');
    xlabel(Ax,'[mm]'); ylabel(Ax,'[mm]');
else %figure was provided in input
    Fig=varargin{1};
    Ax=findobj(Fig,'type','axes');
    cla(Ax);
end

%Initalize Drawing
[w_mm,xiVec_mm]=Compute_w(q,L,100);
hBeam=plot(Ax,xiVec_mm,w_mm,'linewidth',10,'color',[0.4,0,0.4]);

%force
hforce=plot(Ax,[0,0],[0,20],'linewidth',10,'color',[0.5,0.7,1]);

%spring
SpringR=200; SpringN=10;
oSpr=Spring(SpringR,SpringN);
[x,y]=oSpr.getSpr([0,0],[0,-30]);
hSpr=plot(Ax,x,y,'k','linew',2);

%static Drawings
scatter(Ax,L/8*1e3,0,200,[0,0.5,0.5],'filled'); %<------Fix for beam
scatter(Ax,0,-30,200,[0,0.5,0.5],'filled','sq'); %<------Fix for spring

%Initalize text for time
xtext=0.9*Ax.XLim(1)+0.1*Ax.XLim(2);
ytext=0.1*Ax.YLim(1)+0.9*Ax.YLim(2);
hTime=text(Ax,xtext,ytext,'');
%% Storing handles to "figure" and block "UserData"
UserData.Figure = Fig;
UserData.Axes = Ax;
UserData.hBeam=hBeam;
UserData.hforce=hforce;
UserData.hTime = hTime;
UserData.oSpr = oSpr;
UserData.hSpr = hSpr;

%Store in both figure and block
set(gcbh,'UserData',UserData);
end
function [w_mm,xiVec_mm]=Compute_w(q,L,N)
f_psi=@(xi)[1.0;xi.*(2.0./5.0)-1.0;xi.*(-8.0./5.0)+xi.^2.*(8.0./2.5e+1)+1.0;xi.*(1.8e+1./5.0)-xi.^2.*(4.8e+1./2.5e+1)+xi.^3.*(3.2e+1./1.25e+2)-1.0;xi.*(-3.2e+1./5.0)+xi.^2.*(3.2e+1./5.0)-xi.^3.*(2.56e+2./1.25e+2)+xi.^4.*(1.28e+2./6.25e+2)+1.0;xi.*1.0e+1-xi.^2.*1.6e+1+xi.^3.*(2.24e+2./2.5e+1)-xi.^4.*(2.56e+2./1.25e+2)+xi.^5.*1.6384e-1-1.0;xi.*(-7.2e+1./5.0)+xi.^2.*(1.68e+2./5.0)-xi.^3.*2.8672e+1+xi.^4.*1.10592e+1-xi.^5.*1.96608+xi.^6.*1.31072e-1+1.0;xi.*(9.8e+1./5.0)-xi.^2.*6.272e+1+xi.^3.*7.5264e+1-xi.^4.*4.3008e+1+xi.^5.*1.261568e+1-xi.^6.*1.835008+xi.^7.*1.048576e-1-1.0];
xiVec=linspace(0,L,N)'; %col vec
w=zeros(N,1);
for i=1:N
w(i)=f_psi(xiVec(i))'*q;
end
w_mm=w*1e3;
xiVec_mm=xiVec*1e3;
end
%% Unused fcns
function Terminate(block)
end
function Start(block)
Enable=block.InputPort(1).Data(1);
if ~Enable, return, end

end
function CheckPrms(block)
  %can check validity of parameters here
end