classdef NetCDFFile < NetCDFGroup
    % Read and write NetCDF files
    %
    % NetCDF files are a standard file format for reading and writing data.
    % This class is designed to simplify the task of adding new dimensions,
    % variables, and attributes to a NetCDF file compared to using the
    % built-in `ncread` and `ncwrite` functions.
    %
    % Typical usage patterns include:
    %   1) Create/open a file and obtain the root group interface.
    %   2) Define dimensions, then define variables referencing those
    %      dimensions.
    %   3) Attach global and variable attributes, and synchronize/close.
    %
    % ```matlab
    % ncfile = NetCDFFile('myfile.nc')
    %
    % % create two new dimensions and add them to the file
    % x = linspace(0,10,11);
    % y = linspace(-10,0,11);
    % ncfile.addDimension('x',x);
    % ncfile.addDimension('y',y);
    %
    % % Create new multi-dimensional variables, and add those to the file
    % [X,Y] = ncgrid(x,y);
    % ncfile.addVariable(X,{'x','y'});
    % ncfile.addVariable(Y,{'x','y'});
    % ```
    %
    % - Topic: Initializing
    % - Topic: Accessing file properties
    % - Topic: Working with dimensions
    % - Topic: Working with variables
    % - Topic: Working with groups
    % - Topic: Working with global attributes
    %
    % - Declaration: classdef NetCDFFile < NetCDFGroup
    properties
        % File path to the NetCDF dataset.
        %
        % - Topic: Accessing file properties
        path

        % NetCDF format identifier.
        %
        % Default is FORMAT_NETCDF4 for newly created files; for opened files
        % this is overwritten by netcdf.inqFormat.
        % - Topic: Accessing file properties
        format = 'FORMAT_NETCDF4'
    end

    properties (Dependent)
        % File name (base name + extension) derived from path.
        % - Topic: Accessing file properties
        filename
    end
    methods
        function self = NetCDFFile(path,options)
            % Open an existing NetCDF file or create a new file at path.
            %
            % If a file exists at path, it will be opened for reading
            % or read/write (depending on options). If no file exists, the
            % a new file will be created.
            %
            % - Topic: Initializing
            % - Declaration: self = NetCDFFile(path,options)
            % - Parameter path: (string) file path
            % - Parameter options.shouldOverwriteExisting: logical scalar — delete existing file at path before creating a new dataset
            % - Parameter options.shouldUseClassicNetCDF: logical scalar — create classic (non-NetCDF4) file format
            % - Parameter options.shouldReadOnly: logical scalar — open existing file in read-only mode
            % - Returns self: NetCDFFile instance
            arguments
                path (1,1) string {mustBeNonempty}
                options.shouldOverwriteExisting logical = false
                options.shouldUseClassicNetCDF logical = false
                options.shouldReadOnly logical = false
            end

            if isfile(path) && options.shouldOverwriteExisting == 1
                delete(path);
                shouldCreateNewFile = 1;
            elseif ~isfile(path)
                shouldCreateNewFile = 1;
            else
                shouldCreateNewFile = 0;
            end

            if shouldCreateNewFile == 1
                if options.shouldUseClassicNetCDF == 1
                    ncid = netcdf.create(path, bitor(netcdf.getConstant('SHARE'),netcdf.getConstant('WRITE')));
                else  
                    ncid = netcdf.create(path, netcdf.getConstant('NETCDF4'));
                end
            else
                if options.shouldReadOnly
                    ncid = netcdf.open(path);
                else
                    ncid = netcdf.open(path, bitor(netcdf.getConstant('SHARE'),netcdf.getConstant('WRITE')));
                end
            end
            
            self@NetCDFGroup(id=ncid);
            self.path = path;
            self.format = netcdf.inqFormat(self.id);            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Utilities
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function filename = get.filename(self)
            % Return the file name (including extension).
            %
            % - Topic: Accessing file properties
            % - Declaration: filename = get.filename()
            % - Returns filename: text scalar — file name derived from path
            arguments
                self NetCDFFile
            end
            [~,name,ext] = fileparts(self.path);
            filename = name + ext;
        end


        function self = sync(self)
            % Synchronize in-memory changes to disk.
            %
            % Calls netcdf.sync on the underlying file handle.
            %
            % - Topic: Accessing file properties
            % - Declaration: self = sync()
            % - Returns self: NetCDFFile instance
            arguments
                self NetCDFFile
            end
            netcdf.sync(self.id);
        end

        % function self = open(self)
        %     % - Topic: Accessing file properties
        %     self.ncid = netcdf.open(self.path, bitor(netcdf.getConstant('SHARE'),netcdf.getConstant('WRITE')));
        %     self.format = netcdf.inqFormat(self.ncid);
        % end

        function self = close(self)
            % Close the underlying NetCDF file handle.
            %
            % After close, the instance's id is cleared.
            %
            % - Topic: Accessing file properties
            % - Declaration: self = close()
            % - Returns self: NetCDFFile instance
            arguments
                self NetCDFFile
            end
            netcdf.close(self.id);
            self.id = [];
        end

        function delete(self)
            % Destructor: close the file handle if still open.
            %
            % - Topic: Accessing file properties
            % - Declaration: delete()
            arguments
                self NetCDFFile
            end
            if ~isempty(self.id)
                netcdf.close(self.id);
                self.id = [];
                disp('NetCDFFile closed.');
            end
        end

        function ncfile = duplicate(self,path,options)
            % Duplicate the NetCDF dataset to a new file.
            %
            % Creates a new NetCDFFile at the requested path and copies the
            % current group hierarchy, dimensions, variables, and attributes
            % using addDuplicateGroup. Optionally restricts copied data via
            % indexRange.
            %
            % - Topic: Working with groups
            % - Declaration: ncfile = duplicate(path,options)
            % - Parameter path: any — destination file path
            % - Parameter options.shouldOverwriteExisting: logical scalar (1x1) — overwrite destination if it exists
            % - Parameter options.indexRange: dictionary (string -> cell) — per-variable slicing ranges forwarded to addDuplicateGroup
            % - Returns ncfile: NetCDFFile instance — duplicated dataset handle
            arguments
                self NetCDFFile
                path
                options.shouldOverwriteExisting logical = false
                options.indexRange = dictionary(string.empty,cell.empty)
            end
            if options.shouldOverwriteExisting == 1
                if isfile(path)
                    delete(path);
                end
            else
                if isfile(path)
                    error('A file already exists with that name.')
                end
            end
            ncfile = NetCDFFile(path);
            ncfile.addDuplicateGroup(self,shouldAddToSelf=true,indexRange=options.indexRange);
        end
    end
end