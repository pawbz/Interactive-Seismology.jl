"""Read and validate the notebook registries used by every deployment target."""
module DeploymentNotebooks

using YAML
import Pluto

export SOURCE_DIRECTORY, notebook_paths, registry_paths

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

"""Return the static, live, and combined notebook paths, relative to `src/`."""
function notebook_paths(; root::AbstractString=REPOSITORY_ROOT)
    path = joinpath(root, "live-notebooks.yml")
    isfile(path) || error("Notebook registry not found: $(path)")

    registry = YAML.load_file(path)
    sections = get(registry, "sections", nothing)
    sections isa AbstractDict || error("live-notebooks.yml must contain a 'sections:' mapping")

    static = String[]
    live = String[]
    for (section, section_entries) in sections
        section_entries isa AbstractDict || begin
            section_entries isa AbstractVector || error("live-notebooks.yml: section $(repr(section)) must be either a list or an object with 'live'/'static' fields")
            for entry in section_entries
                entry isa AbstractString || error("live-notebooks.yml: notebook entries must be strings")
                normalized = normalize_notebook_entry(entry, "live-notebooks.yml"; root)
                push!(live, normalized)
            end
            continue
        end

        for kind in ("live", "static")
            entries = get(section_entries, kind, nothing)
            entries === nothing && continue
            entries isa AbstractVector || error("live-notebooks.yml: section $(repr(section)) field $(repr(kind)) must contain a list")
            for entry in entries
                entry isa AbstractString || error("live-notebooks.yml: notebook entries must be strings")
                normalized = normalize_notebook_entry(entry, "live-notebooks.yml"; root)
                if kind == "live"
                    push!(live, normalized)
                else
                    push!(static, normalized)
                end
            end
        end
    end

    length(unique(live)) == length(live) || error("live-notebooks.yml contains a live notebook more than once")
    length(unique(static)) == length(static) || error("live-notebooks.yml contains a static notebook more than once")
    overlap = intersect(static, live)
    isempty(overlap) || error("A notebook must be static or live, not both: $(join(overlap, ", "))")
    (; static, live, all=vcat(static, live))
end

end
