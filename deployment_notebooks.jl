"""Read and validate the notebook registries used by every deployment target."""
module DeploymentNotebooks

using YAML
import Pluto

export SOURCE_DIRECTORY, notebook_entries, notebook_paths, registry_paths

const REPOSITORY_ROOT = @__DIR__
const SOURCE_DIRECTORY = joinpath(REPOSITORY_ROOT, "src")

function normalize_notebook_entry(entry::AbstractString, filename::AbstractString; root::AbstractString=REPOSITORY_ROOT)
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
    relative_path
end

function registry_paths(filename::AbstractString; root::AbstractString=REPOSITORY_ROOT)
    path = joinpath(root, filename)
    isfile(path) || error("Notebook registry not found: $(path)")

    registry = YAML.load_file(path)
    sections = get(registry, "sections", nothing)
    sections isa AbstractDict || error("$(filename) must contain a 'sections:' mapping")

    paths = String[]
    for (section, section_entries) in sections
        if section_entries isa AbstractDict
            for kind in ("live", "static")
                entries = get(section_entries, kind, nothing)
                entries === nothing && continue
                entries isa AbstractVector || error("$(filename): section $(repr(section)) field $(repr(kind)) must contain a list")
                for entry in entries
                    entry isa AbstractString || error("$(filename): notebook entries must be strings")
                    push!(paths, normalize_notebook_entry(entry, filename; root))
                end
            end
        elseif section_entries isa AbstractVector
            for entry in section_entries
                entry isa AbstractString || error("$(filename): notebook entries must be strings")
                push!(paths, normalize_notebook_entry(entry, filename; root))
            end
        else
            error("$(filename): section $(repr(section)) must be either a list or an object with 'live'/'static' fields")
        end
    end

    length(unique(paths)) == length(paths) || error("$(filename) contains a notebook more than once")
    paths
end

"""Return the unified registry entries, preserving the order within each YAML list."""
function notebook_entries(; root::AbstractString=REPOSITORY_ROOT)
    filename = "live-notebooks.yml"
    path = joinpath(root, filename)
    isfile(path) || error("Notebook registry not found: $(path)")

    registry = YAML.load_file(path)
    sections = get(registry, "sections", nothing)
    sections isa AbstractDict || error("$(filename) must contain a 'sections:' mapping")

    entries = NamedTuple{(:section, :kind, :path), Tuple{String, Symbol, String}}[]
    for (section, section_entries) in sections
        if section_entries isa AbstractDict
            for kind in ("live", "static")
                listed = get(section_entries, kind, nothing)
                listed === nothing && continue
                listed isa AbstractVector || error("$(filename): section $(repr(section)) field $(repr(kind)) must contain a list")
                for entry in listed
                    entry isa AbstractString || error("$(filename): notebook entries must be strings")
                    push!(entries, (section=String(section), kind=Symbol(kind), path=normalize_notebook_entry(entry, filename; root)))
                end
            end
        elseif section_entries isa AbstractVector
            for entry in section_entries
                entry isa AbstractString || error("$(filename): notebook entries must be strings")
                push!(entries, (section=String(section), kind=:live, path=normalize_notebook_entry(entry, filename; root)))
            end
        else
            error("$(filename): section $(repr(section)) must be either a list or an object with 'live'/'static' fields")
        end
    end

    paths = getproperty.(entries, :path)
    length(unique(paths)) == length(paths) || error("$(filename) contains a notebook more than once")
    entries
end

"""Return the static, live, and combined notebook paths, relative to `src/`."""
function notebook_paths(; root::AbstractString=REPOSITORY_ROOT)
    entries = notebook_entries(; root)
    static = [entry.path for entry in entries if entry.kind == :static]
    live = [entry.path for entry in entries if entry.kind == :live]
    (; static, live, all=vcat(static, live))
end

end
