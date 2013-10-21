
function init = struct_override(init, init_)
    if isstruct(init_)
        for fld = fieldnames(init_)'
            init.(fld{1}) = init_.(fld{1});
        end
    end   
end