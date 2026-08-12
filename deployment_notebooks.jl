"""Read and validate the notebook registries used by every deployment target."""
module DeploymentNotebooks

using YAML
import Pluto

export SOURCE_DIRECTORY, notebook_paths, registry_paths

const REPOSITORY_ROOT = @__DIR__
const SOURCE_DIRECTORY = joinpath(REPOSITORY_ROOT, "src")

function registry_paths(filename::AbstractString; root::AbstractString=REPOSITORY_ROOT)
    path = joinpath(root, filename)
    isfile(path) || error("Notebook registry not found: $(path)")

    registry = YAML.load_file(path)
    sections = get(registry, "sections", nothing)
    sections isa AbstractDict || error("$(filename) must contain a 'sections:' mapping")

    paths = String[]
    for (section, entries) in sections
        entries isa AbstractVector || error("$(filename): section $(repr(section)) must contain a list")
        for entry in entries
            entry isa AbstractString || error("$(filename): notebook entries must be strings")
            relative_path = normpath(String(entry))
            isabspath(relative_path) && error("$(filename): absolute paths are not allowed: $(entry)")
            (relative_path == ".." || startswith(relative_path, ".." * Base.Filesystem.path_separator)) &&
                error("$(filename): paths must stay inside src: $(entry)")
            endswith(relative_path, ".jl") || error("$(filename): notebook must end in .jl: $(entry)")
            notebook_path = joinpath(root, "src", relative_path)
            isfile(notebook_path) ||
                error("$(filename): notebook does not exist: src/$(relative_path)")
            Pluto.is_pluto_notebook(notebook_path) ||
                error("$(filename): file is not a Pluto notebook: src/$(relative_path)")
            push!(paths, relative_path)
        end
    end

    length(unique(paths)) == length(paths) || error("$(filename) contains a notebook more than once")
    paths
end

"""Return the static, live, and combined notebook paths, relative to `src/`."""
function notebook_paths(; root::AbstractString=REPOSITORY_ROOT)
    static = registry_paths("static-notebooks.yml"; root)
    live = registry_paths("live-notebooks.yml"; root)
    overlap = intersect(static, live)
    isempty(overlap) || error("A notebook must be static or live, not both: $(join(overlap, ", "))")
    (; static, live, all=vcat(static, live))
end

end
