classdef FigureSaver

    properties
        folder = 'figures';
        % units = 'inches';
        % paper_width = 11;
        % paper_height = 8.5;
        % padding = 0.3;
    end

    methods
        function obj = FigureSaver(folder)
            arguments
                folder = 'figures/'
            end

            obj.folder = folder;
        end

        function SavePDF(obj, figure_name)
            % SavePDF(figure_name)
            % where figure_name is a string without the .pdf extension

            figure_path = fullfile(obj.folder, figure_name + ".pdf");

            print('-dpdf', '-bestfit', figure_path)

        end

        function SavePDF_fill(obj, figure_name)

            figure_path = fullfile(obj.folder, figure_name + ".pdf");

            print('-dpdf', '-fillpage', figure_path)


        end

    end
end


% exportgraphics(fig, figure_path, ...
%     'units' , obj.units, ...
%     'width' , obj.paper_width, ...
%     'height' , obj.paper_height,...
%     'padding' , obj.padding, ...
%     'contenttype', 'vector');
