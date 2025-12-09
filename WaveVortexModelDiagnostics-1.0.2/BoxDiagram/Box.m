%% File: Box.m
classdef Box < handle
    %BOX Encapsulates a labeled rectangle that can be drawn on an axes.
    properties
        Label string                         % text inside the box
        Position (1,2) double % lower-left corner [x y]
        Size (1,2) double % [width height]
        Sublabel  string   = ""
        CornerRadius double {mustBeNonnegative} = 0 % 0 = square corners
        FaceColor (1,3) double = [1 1 1]
        EdgeColor (1,3) double = [0 0 0]
        FontSize  double {mustBePositive} = 16
        FontSizeSublabel  double {mustBePositive} = 12
    end

    methods
        function obj = Box(label,pos,size,options)
            %BOX Construct a new Box
            arguments
                label string
                pos  (1,2) double = [nan nan]
                size (1,2) double = [nan nan]
                options.Sublabel  string   = ""
                options.CornerRadius double = 0
                options.FaceColor   (1,3) double = [1 1 1]
                options.EdgeColor   (1,3) double = [0 0 0]
                options.FontSize    double = 16
                options.FontSizeSublabel    double = 12
            end
            obj.Label        = label;
            obj.Position     = pos;
            obj.Size         = size;
            obj.CornerRadius = options.CornerRadius;
            obj.FaceColor    = options.FaceColor;
            obj.EdgeColor    = options.EdgeColor;
            obj.FontSize     = options.FontSize;
            obj.Sublabel     = options.Sublabel;
            obj.FontSizeSublabel = options.FontSizeSublabel;
        end
        function setPosition(obj,pos)
            %SETPOSITION Assign [x y] lower-left corner after construction.
            obj.Position = pos;
        end
        
        function shp = polyshape(self)
            xmin = self.Position(1);
            xmax = self.Position(1) + self.Size(1);
            ymin = self.Position(2);
            ymax = self.Position(2) + self.Size(2);
            shp = polyshape([xmin xmax xmax xmin], [ymin ymin ymax ymax]);
        end

        function h = draw(obj,ax)
            %DRAW Render the box at its current Position.
            if isempty(obj.Position)
                error("Box:PositionUnset","Position must be set before drawing.");
            end
            if nargin<2 || isempty(ax), ax = gca; end
            x = obj.Position(1); y = obj.Position(2);
            w = obj.Size(1);     hgt = obj.Size(2);
            if obj.CornerRadius>0
                h = rectangle(ax,'Position',[x y w hgt], ...
                    'Curvature',[obj.CornerRadius obj.CornerRadius*w/hgt], ...
                    'FaceColor',obj.FaceColor,'EdgeColor',obj.EdgeColor,'LineWidth',1.2);
            else
                h = rectangle(ax,'Position',[x y w hgt], ...
                    'FaceColor',obj.FaceColor,'EdgeColor',obj.EdgeColor,'LineWidth',1.2);
            end
            xc = x+w/2; yc = y+hgt/2;
            if strlength(obj.Sublabel)>0
                hMain = text(ax,xc,yc,obj.Label, HorizontalAlignment='center',VerticalAlignment='bottom',FontSize=obj.FontSize,Interpreter='none',Visible='off');
                extMain = get(hMain,"Extent");  % [x y width height]
                hMainHeight = extMain(4);
                % next, draw sublabel invisibly (if any) to get its extent
                if strlength(obj.Sublabel)>0
                    hSub = text(ax, xc, yc, obj.Sublabel,HorizontalAlignment="center", VerticalAlignment="top",FontSize=obj.FontSizeSublabel,Visible="off");
                    extSub   = get(hSub,"Extent");
                    hSubHeight = extSub(4);
                else
                    hSubHeight = 0;
                end

                % compute common y so that block's center = yc
                yCommon = yc + (hSubHeight - hMainHeight)/2;

                % now make them visible at the correct spot
                set(hMain, Position=[xc yCommon 0],VerticalAlignment="bottom",Visible="on");

                if strlength(obj.Sublabel)>0
                    set(hSub, Position=[xc yCommon 0],VerticalAlignment="top",Visible="on");
                end
            else
                text(ax,xc,yc,obj.Label, HorizontalAlignment='center',VerticalAlignment='middle',FontSize=obj.FontSize,Interpreter='none');
            end
        end

        function [xc,yc,hw,hh] = getGeometry(obj)
            % Return centre and half-sizes (for arrow maths).
            if isempty(obj.Position)
                error("Box:PositionUnset","Position must be set before geometry queries.");
            end
            hw = obj.Size(1)/2; hh = obj.Size(2)/2;
            xc = obj.Position(1) + hw;
            yc = obj.Position(2) + hh;
        end

        function [xe,ye] = edgeIntersect(obj,theta)
            %EDGEINTERSECT Intersection of ray angle theta from centre with edge.
            [xc,yc,hw,hh] = obj.getGeometry();
            dx = cos(theta); dy = sin(theta);
            if abs(dx)<1e-12, dx = 1e-12*sign(dx+eps); end
            if abs(dy)<1e-12, dy = 1e-12*sign(dy+eps); end
            tx = hw/abs(dx); ty = hh/abs(dy);
            t  = min(tx,ty);
            xe = xc + t*dx;
            ye = yc + t*dy;
        end
    end
end