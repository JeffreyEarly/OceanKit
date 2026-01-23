classdef BoxDiagram < handle
    properties
        rows cell
        arrows Arrow
        arrowsToDraw cell = {}
        RowLabels
        Title string = ""
        TitleFontSize = 16;
        BackgroundColor (1,3) double = [1 1 1]
    end

    methods

        function self = BoxDiagram(row1,row2,row3, arrows, options)
            arguments
                row1 Box
                row2 Box
                row3 Box
                arrows Arrow
                options.nRow2Subrows = 2
                options.RowLabels (1,3) string = ["Sources" "Reservoirs" "Sinks"]
                options.RowSublabels (1,3) string = strings(1,3)
                options.BoxSize (1,2) double {mustBePositive} = [2 1]
                options.Dx double {mustBeNonnegative} = 1
                options.Dy double {mustBeNonnegative} = 2
                options.Origin (1,2) double = [0 9]
                options.Title string = ""
                options.BackgroundColor (1,3) double = [1 1 1]
                options.BaseLineWidth double {mustBePositive} = 1.5
            end

            for i=1:length(options.RowLabels)
                self.RowLabels(i).label = options.RowLabels(i);
                self.RowLabels(i).labelFontSize = 16;
                self.RowLabels(i).sublabel = options.RowSublabels(i);
                self.RowLabels(i).sublabelFontSize = 12;
            end

            self.Title = options.Title;
            self.BackgroundColor = options.BackgroundColor;

            self.rows{1} = row1;
            boxesPerRow = ceil(length(row2)/options.nRow2Subrows);
            iRow = 2;
            for i=1:options.nRow2Subrows
                if i==options.nRow2Subrows
                    indices = ((i-1)*boxesPerRow+1):length(row2);
                else
                    indices = ((i-1)*boxesPerRow+1):i*boxesPerRow;
                end
                self.rows{iRow} = row2(indices); iRow = iRow+1;
            end
            self.rows{iRow} = row3;

            self.arrows = arrows;

            self.layoutDiagram(BoxSize=options.BoxSize,Dx=options.Dx,Dy=options.Dy,Origin=options.Origin,BaseLineWidth=options.BaseLineWidth)
        end

        function layoutDiagram(self,options)
            arguments
                self BoxDiagram
                options.BoxSize (1,2) double {mustBePositive} = [2 1]
                options.Dx double {mustBeNonnegative} = 1
                options.Dy double {mustBeNonnegative} = 2
                options.Origin (1,2) double = [0 9]
                options.BaseLineWidth double {mustBePositive} = 1.5
                options.MaxLineWidth double {mustBePositive} = 5
            end

            self.layoutBoxes(BoxSize=options.BoxSize,Dx=options.Dx,Dy=options.Dy,Origin=options.Origin);
            
            self.layoutArrows(BaseLineWidth=options.BaseLineWidth,MaxLineWidth=options.MaxLineWidth);

            lblOffset = 2.0;
            for r = 1:3
                if r==1
                    yMid  = self.rows{r}(1).Position(2) + self.rows{r}(1).Size(2)/2;
                elseif r==2
                    yMid  = self.rows{r+1}(1).Position(2) + self.rows{r+1}(1).Size(2)*3/2;
                elseif r==3
                    yMid  = self.rows{r+1}(1).Position(2) + self.rows{r+1}(1).Size(2)/2;
                end
                xLeft = options.Origin(1) - lblOffset;

                self.RowLabels(r).xLeft = xLeft;
                self.RowLabels(r).yMid = yMid;
            end
        end

        function layoutBoxes(self,options)
            arguments
                self BoxDiagram
                options.BoxSize (1,2) double {mustBePositive} = [2 1]
                options.Dx double {mustBeNonnegative} = 1
                options.Dy double {mustBeNonnegative} = 2
                options.Origin (1,2) double = [0 9]
            end

            nRows = length(self.rows);
            w      = options.BoxSize(1); h = options.BoxSize(2);
            dx     = options.Dx;        dy = options.Dy;
            x0     = options.Origin(1); y0 = options.Origin(2);

            maxN   = max(cellfun(@numel,self.rows));

            layout_style = "fixed-spacing-centered";
            layout_style = "max-spacing";
            if layout_style == "fixed-spacing-centered"
                % Assign positions to boxes
                for r = 1:nRows
                    thisRow = self.rows{r}; n = numel(thisRow);
                    xs = x0 + ((maxN - n)*(w + dx))/2 + (0:n-1)*(w + dx);
                    y  = y0 - (r-1)*(h + dy);
                    for k = 1:n
                        thisRow(k).Size = [w h];
                        thisRow(k).setPosition([xs(k) y]);
                    end
                end
            elseif layout_style == "max-spacing"
                maxWidth = (maxN-1)*(w + dx) + w;
                for r = 1:nRows
                    thisRow = self.rows{r}; n = numel(thisRow);
                    if n > 1
                        dx = (maxWidth - w)/(n-1)-w;
                        xs = x0 + (0:n-1)*(w + dx);
                    else
                        xs = x0 + (maxWidth-w)/2;
                    end
                    y  = y0 - (r-1)*(h + dy);
                    for k = 1:n
                        thisRow(k).Size = [w h];
                        thisRow(k).setPosition([xs(k) y]);
                    end
                end
            end
        end

        function layoutArrows(self,options)
            arguments
                self BoxDiagram
                options.BaseLineWidth double {mustBePositive} = 1.5
                options.MaxLineWidth double {mustBePositive} = 5
            end

            % Determine line‑width scale from Magnitude field
            mags   = [self.arrows.Magnitude];
            maxMag = max(mags);
            % minMag = max(min(mags),0.01);
            boxes = [self.rows{:}];
            for i = 1:length(self.arrows)
                a = self.arrows(i);
                a.LineWidth = options.BaseLineWidth + (options.MaxLineWidth-options.BaseLineWidth)*(a.Magnitude/maxMag);
                if BoxDiagram.arrowIntersectsAnyBox(a, boxes, find(ismember(boxes,[a.Source a.Target])))
                    disp("Arrow " + string(i) + " with label " + a.Label + " intersects a box.");
                end
            end
        end

        function fig = draw(self,options)
            arguments
                self BoxDiagram
                options.visible = "on"
            end

            fig = figure('Color',self.BackgroundColor,'Position',[100 100 800 600],Visible=options.visible);
            ax   = axes('Parent',fig); hold(ax,'on');  axis(ax,'off'); axis(ax,'equal');

            % Draw boxes
            nRows = length(self.rows);
            for r = 1:nRows, arrayfun(@(b) b.draw(ax), self.rows{r}); end

            % Draw arrows
            arrayfun(@(a) a.draw(ax), self.arrows);
            % cellfun(@(a) a.draw(ax), self.arrowsToDraw);

            for r = 1:3
                text(ax,self.RowLabels(r).xLeft,self.RowLabels(r).yMid,self.RowLabels(r).label,'FontSize',self.RowLabels(r).labelFontSize,'FontWeight','bold', ...
                    'Rotation',90,'HorizontalAlignment','center','VerticalAlignment','middle');
                if strlength(self.RowLabels(r).sublabel)
                    text(ax,self.RowLabels(r).xLeft + 0.5,self.RowLabels(r).yMid,self.RowLabels(r).sublabel,'FontSize',self.RowLabels(r).sublabelFontSize,'FontWeight','normal', ...
                        'Rotation',90,'HorizontalAlignment','center','VerticalAlignment','middle');
                end
            end

            if strlength(self.Title)>0
                title(ax,self.Title,'FontSize',self.TitleFontSize,'FontWeight','bold', 'Interpreter', 'none');
            end

            ax.Clipping = "off";
        end
    end

    methods (Static)
        function tf = arrowIntersectsAnyBox(arrow, boxes, ignoreIdx)
            %ARROWINTERSECTSANYBOX true if segment p0->p1 intersects any box rect.
            % p0, p1 : 1x2 numeric [x y]
            % boxes  : array of Box objects
            % ignoreIdx : indices to skip (e.g., source/target boxes)

            if nargin < 3, ignoreIdx = []; end
            arrowLine = arrow.polyshape;
            n = numel(boxes);
            for k = 1:n
                if any(k == ignoreIdx)
                    continue;
                end
                if area(intersect(boxes(k).polyshape,arrowLine)) > 0
                    tf = true; return;
                end
            end

            tf = false;
        end
    end

end