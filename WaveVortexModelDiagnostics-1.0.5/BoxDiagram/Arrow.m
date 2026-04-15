classdef Arrow < handle
    %ARROW Connect two Box objects.  Magnitude is a data value which the
    %layout engine can translate into a line‑width; LineWidth is filled in
    %later.  Arrowhead geometry matches the original ShowBoxPlot script.

    properties
        Source Box
        Target Box
        Label string = ""
        Magnitude double {mustBeNonnegative} = 1   % data value
        LineWidth double = NaN                     % data-space width, set later
        Color (1,3) double = [0 0 0]
        LabelOffset double = 0.5
        FontSize double {mustBePositive} = 10
        intermediatePoints = []
        LabelPosition = []
        LabelAlpha = 1.0;
    end

    properties (Constant, Access=private)
        BaseWidthPt = 6;    % head size in printer points when LineWidth = 1 pt
        WidthRatio  = 1.1;  % head base width : head length
    end

    methods
        function obj = Arrow(src,dst,options)
            arguments
                src Box
                dst Box
                options.Label string = ""
                options.Magnitude double = 1
                options.LineWidth double = NaN
                options.Color (1,3) double = [0 0 0]
                options.LabelOffset double = 0.5
                options.FontSize double = 10
            end
            obj.Source      = src;
            obj.Target      = dst;
            obj.Label       = options.Label;
            obj.Magnitude   = options.Magnitude;
            obj.LineWidth   = options.LineWidth;
            obj.Color       = options.Color;
            obj.LabelOffset = options.LabelOffset;
            obj.FontSize    = options.FontSize;
        end

        function pt = sourcePoint(obj)
            [xc1,yc1] = obj.Source.getGeometry();
            pt = [xc1,yc1];
        end

        function pt = targetPoint(obj)
            [xc2,yc2] = obj.Target.getGeometry();
            pt = [xc2,yc2];
        end

        function shp = polyshape(self)
            shp = polybuffer([self.sourcePoint; self.targetPoint],'lines',0.1);
        end

        function draw(self,ax)
            if isnan(self.LineWidth) || self.LineWidth<=0
                error("Arrow:LineWidthUnset","LineWidth must be assigned before drawing");
            end
            if nargin<2 || isempty(ax), ax = gca; end

            points = cat(1,self.sourcePoint,self.intermediatePoints,self.targetPoint);
            i=0;
            for j=1:(size(points,1)-2)
                xc1=points(j,1); yc1 = points(j,2);
                x2=points(j+1,1); y2 = points(j+1,2);
                ang = atan2(y2-yc1, x2-xc1);
                if j==1
                    [x1,y1] = self.Source.edgeIntersect(ang);
                else
                    x1 = xc1; y1 = yc1;
                end
                plot(ax,[x1 x2],[y1 y2],'Color',self.Color,'LineWidth',self.LineWidth);
                i=j;
            end

            %  Geometry
            xc1=points(i+1,1); yc1 = points(i+1,2);
            xc2=points(i+2,1); yc2 = points(i+2,2);
            ang = atan2(yc2-yc1, xc2-xc1);
            if i+1==1
                [x1,y1] = self.Source.edgeIntersect(ang);
            else
                x1 = xc1; y1 = yc1;
            end
            [x2,y2] = self.Target.edgeIntersect(ang+pi);

            %  Head dimensions in data units (proportional to LineWidth)
            headLen = (Arrow.BaseWidthPt/72) * self.LineWidth; % approx
            headWidth = Arrow.WidthRatio * headLen;

            %  Shaft (trimmed before head)
            x2b = x2 - headLen*cos(ang);
            y2b = y2 - headLen*sin(ang);
            plot(ax,[x1 x2b],[y1 y2b],'Color',self.Color,'LineWidth',self.LineWidth);

            %  Head triangle
            left  = [x2b + headWidth/2*sin(ang), y2b - headWidth/2*cos(ang)];
            right = [x2b - headWidth/2*sin(ang), y2b + headWidth/2*cos(ang)];
            patch(ax,[x2 left(1) right(1)],[y2 left(2) right(2)],self.Color,'EdgeColor',self.Color);

            %  Optional label
            if strlength(self.Label)>0
                if isempty(self.LabelPosition)
                tx = x1 + self.LabelOffset*(x2 - x1);
                ty = y1 + self.LabelOffset*(y2 - y1) + 0.025;
                else
                    tx = self.LabelPosition(1); ty = self.LabelPosition(2);
                end
                text(ax,tx,ty,self.Label,'HorizontalAlignment','center', ...
                    'VerticalAlignment','middle','FontSize',self.FontSize, ...
                    'BackgroundColor',[1 1 1 self.LabelAlpha],'Margin',0.5,'Interpreter','none');
            end
        end
    end
end