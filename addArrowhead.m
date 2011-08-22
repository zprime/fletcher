% Add an arrowhead to a line at a point
%
% function addArrowhead(hline,points,OPTIONS)
%
% Add an arrowhead to a line. The arrowhead remains a constant size and
% aspect ratio independent of the data size and aspect ratio.  The
% arrowhead also curves with the data line, making arrows on curved lines
% look nice.
%
% INPUTS:
%   hline (required)     - handle to the line to put the arrowhead on.
%   points (required)    - point on the line. Can either be an array
%                          indices or the string 'end'.
%   OPTIONS: (optional key-value pairs)
%   Key         | Value
%   -----------------------------------------------------------------------
%   'direction' | ['start'|'end'] point to the start or the end. By default
%               | if point=1, direction='start', otherwise direction='end'.
%               |
%   'angle'     | Arrowhead half-angle in degrees. The default is 10 deg.
%               |
%   'length'    | Arrowhead length relative the the linewidth of the parent
%               | line. For example if 'length'=20 and hline linewidth=0.5,
%               | then the arrowhead will be 10pts long. The default if 10.
%               |
%   'style'     | ['filled','lines'] Arrowhead style, either a 'filled'
%               | block or two lines. The default is 'block'
%               |
%   'density'   | Number of points used to contruct each edge of the arrow.
%               | The default is 10.
%
% EXAMPLE:
%   x = linspace(1,10,101);
%   hl = plot(x,sin(2*x));
%   addArrowhead(hl,1,'style','lines','angle',30);
%   addArrowhead(hl,'end','sytle','lines','angle',30);
%   addArrowhead(hl,11:10:91);
%
% 20-Aug-2011 v0.1dev
%
% License and Copyright statement appended to the source

function addArrowhead(hline,points,varargin)

% Parse the inputs.
ip = inputParser;
ip.addRequired('hline', @(x) isprop(x,'type') && strcmpi( get(x,'type'),'line' ) );
ip.addRequired('points', @(x) (isnumeric(x) && (all(x>0) && all(x<=length( get(hline,'xdata'))))) ...
  || (strcmpi(x,'end')) );
ip.addParamValue('direction', 'end', @(x) any( strcmpi(x,{'start','end'}) ) );
ip.addParamValue('angle', 10, @(x) isnumeric(x) && (length(x)==1) && x>0 && x<180 );
ip.addParamValue('length', 10, @(x) isnumeric(x) && (length(x)==1) && (x>0) );
ip.addParamValue('style', 'filled', @(x) any(strcmpi(x, {'filled', 'lines'} )) );
ip.addParamValue('density',10, @(x) isnumeric(x) && (length(x)==1) && (x>0) );

ip.parse( hline, points, varargin{:} );

% If point==1, use the direction 'start'
arrowdata.dir = ip.Results.direction;
if any( strcmpi( ip.UsingDefaults, 'direction' ) )
  if all(ip.Results.points == 1)
    arrowdata.dir = 'start';
  end
end

% Points (numerical values)
if strcmp( ip.Results.points ,'end')
  points = length( get(hline,'xdata') );
else
  points = ip.Results.points;
end

% Angle
arrowdata.ang = ip.Results.angle;

% Relative length
arrowdata.rlen = ip.Results.length;

% Arrow density
arrowdata.dens = ip.Results.density;

% Create a hggroup if necessary
hlp = get(hline,'parent');
if strcmpi( get(hlp,'type'), 'hggroup' )
  hg = hlp;
  ha = get(hg,'parent');
else
  hg = hggroup('parent',hlp);
  ha = hlp;
  set(hline,'parent',hg);
end

% Store the axes and line handles in the arrow
arrowdata.hl = hline;
arrowdata.ha = ha;

% Are the arrowheads filled or lines?
arrowdata.filled = strcmpi( ip.Results.style, 'filled' );

% x and y data from the line
arrowdata.xd = get(hline,'xdata');
arrowdata.yd = get(hline,'ydata');

% Loop through and add the arrowheads
for ii=1:length(points)
  % Data index
  arrowdata.n = points(ii);
  
  % X and Y data differences
  if strcmpi(arrowdata.dir,'start')
    arrowdata.dx = circshift( arrowdata.xd, [1 -1] ) - arrowdata.xd;
    arrowdata.dy = circshift( arrowdata.yd, [1 -1] ) - arrowdata.yd;
  else
    arrowdata.dx = arrowdata.xd - circshift( arrowdata.xd, [1 1] );
    arrowdata.dy = arrowdata.yd - circshift( arrowdata.yd, [1 1] );
  end
  %xds = circshift( arrowdata.xd, [1,-points(ii)+1] );
  %yds = circshift( arrowdata.yd, [1,-points(ii)+1] );
  %arrowdata.dx = circshift( [0, diff(xds)], [1,points(ii)-1] );
  %arrowdata.dy = circshift( [0, diff(yds)], [1,points(ii)-1] );
  
  % Plot the base point, then call the redraw subfunction later
  ac = get(hline,'color');
  anp = get(ha,'nextplot');
  set(ha,'nextplot','add');
  if arrowdata.filled
    harr = fill( arrowdata.xd(points(ii)), arrowdata.yd(points(ii)), ...
       ac, 'edgecolor', ac, 'parent', hg );
  else
    harr = plot( arrowdata.xd(points(ii)), arrowdata.yd(points(ii)), ...
      'color', ac, 'parent', hg );
  end
  set(ha,'nextplot',anp);
  
  % Create the listeners
  arrowdata.lxl = handle.listener(handle(ha),findprop(handle(ha),'Xlim'),'PropertyPostSet',{@arrowcallback,harr});
  arrowdata.lyl = handle.listener(handle(ha),findprop(handle(ha),'Ylim'),'PropertyPostSet',{@arrowcallback,harr});
  arrowdata.llw = handle.listener(handle(hline),findprop(handle(hline),'LineWidth'),'PropertyPostSet',{@arrowcallback,harr});
  arrowdata.lc = handle.listener(handle(hline),findprop(handle(hline),'Color'),'PropertyPostSet',{@arrowcallback,harr});
  
  % Now store all the data in the arrow handle
  setappdata( harr, 'ArrowData', arrowdata );
  try
    redrawarrow( harr );
  catch err
    fprintf(1,'The following error occured while drawing an arrow:\n');
    delete( harr );
    rethrow(err);
  end
end % of for ii=1:length(points)

end % of addArrowhead.m

function arrowcallback( eventSrc, eventData, harr )
  arrowdata = getappdata( harr, 'ArrowData' );
  switch lower( get(eventSrc,'name') )
    case {'xlim','ylim','linewidth'}
      redrawarrow( harr );
      
    % If it is color that has changed, simply update that.
    case 'color'
      lc = get(arrowdata.hl,'color');
      if arrowdata.filled
        set( harr, 'edgecolor', lc, 'facecolor', lc );
      else
        set( harr, 'color', lc );
      end
      
    otherwise
      fprintf(1,'Unexpected callback\n');
  end
end % of arrowcallback

% Do the actual arrow redrawing
function redrawarrow( harr )
  arrowdata = getappdata( harr, 'ArrowData' );
  % Get the axes position in points
  hau = get( arrowdata.ha, 'units' );
  set( arrowdata.ha, 'units', 'points' );
  appt = get( arrowdata.ha, 'position' );
  set( arrowdata.ha, 'units', hau );
      
  % Now calculate the x and y scales
  xlp = get( arrowdata.ha, 'xlim' );
  ylp = get( arrowdata.ha, 'ylim' );
  xs = appt(3) / ( xlp(2) - xlp(1) );
  ys = appt(4) / ( ylp(2) - ylp(1) );
  
  % Pre-scale the data
  xd = xs*arrowdata.xd;
  yd = ys*arrowdata.yd;
  
  % Positive data rotation
  xp = ( xd - xd(arrowdata.n)  ).*cosd( arrowdata.ang ) ...
    - ( yd - yd(arrowdata.n) ).*sind( arrowdata.ang ) ...
    + xd(arrowdata.n);
  yp = ( xd - xd(arrowdata.n) ).*sind( arrowdata.ang ) ...
    + ( yd - yd(arrowdata.n) ).*cosd( arrowdata.ang ) ...
    + yd(arrowdata.n);
  
  % Negative data rotation
  xn = ( xd - xd(arrowdata.n)  ).*cosd( arrowdata.ang ) ...
    + ( yd - yd(arrowdata.n) ).*sind( arrowdata.ang ) ...
    + xd(arrowdata.n);
  yn = - ( xd - xd(arrowdata.n) ).*sind( arrowdata.ang ) ...
    + ( yd - yd(arrowdata.n) ).*cosd( arrowdata.ang ) ...
    + yd(arrowdata.n);
  
  % Calculate curve distances in pts
  distd = sqrt( (xs.*arrowdata.dx).^2 + (ys.*arrowdata.dy).^2 );
  ni = linspace(0, arrowdata.rlen*get(arrowdata.hl,'linewidth') ,arrowdata.dens);
  if strcmpi( arrowdata.dir, 'end' )
    dists = cumsum( distd(arrowdata.n:-1:1) );
    nd = arrowdata.n - interp1( [0,dists], 0:length(dists), ni, 'linear', 'extrap' );
  else
    dists = cumsum( distd(arrowdata.n:1:end) );
    nd = arrowdata.n + interp1( [0,dists], 0:length(dists), ni, 'linear', 'extrap' );
  end
  
  % X and Y data interpolation
  ln = 1:length(xp);
  x_pos = interp1( ln, xp, nd, 'linear', 'extrap' );
  x_neg = interp1( ln, xn, nd, 'linear', 'extrap' );
  y_pos = interp1( ln, yp, nd, 'linear', 'extrap' );
  y_neg = interp1( ln, yn, nd, 'linear', 'extrap' );
  
  % Update the arrow data
  if arrowdata.filled
    % Find the mid-point to start at (to make the arrow corners nice)
    x_mid = (x_pos(end)-x_neg(end))/2 + x_neg(end);
    y_mid = (y_pos(end)-y_neg(end))/2 + y_neg(end);
    set( harr, 'xdata', [x_mid,x_pos(end:-1:1),x_neg,x_mid]./xs, 'ydata', [y_mid,y_pos(end:-1:1),y_neg,y_mid]./ys );
  else
    set( harr, 'linewidth', get( arrowdata.hl, 'linewidth' ) );
    set( harr, 'xdata', [x_pos(end:-1:1),x_neg]./xs, 'ydata', [y_pos(end:-1:1),y_neg]./ys );
  end
  
end % of redrawarrow

% Copyright (c) 2011, Zebb Prime and The University of Adelaide
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of the organization nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL ZEBB PRIME OR THE UNIVERSITY OF ADELAIDE BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
% AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.