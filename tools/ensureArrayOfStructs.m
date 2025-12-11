function s = ensureArrayOfStructs(s, fieldNames)
%ENSUREARRAYOFSTRUCTS Force certain fields to encode as JSON arrays of objects.
%
% jsondecode turns a JSON array of objects into a struct array. For a
% single-element array, that becomes a scalar struct, and jsonencode would
% normally output a bare object instead of [ { ... } ].
%
% Converting a struct (or struct array) to a cell array of structs forces
% jsonencode to emit a JSON array, even for one element.

    for f = fieldNames(:).'
        name = char(f);
        if ~isfield(s, name)
            continue
        end

        val = s.(name);

        if isempty(val)
            % Encode as [] rather than {} or nothing
            s.(name) = {};
        elseif isstruct(val)
            % struct or struct array -> cell array of structs
            s.(name) = num2cell(val);
        elseif iscell(val)
            % Assume caller already put it in a cell array; do nothing
        end
    end
end