"""Read and validate the notebook registries used by every deployment target."""
module DeploymentNotebooks

using YAML
import Pluto

export SOURCE_DIRECTORY, notebook_entries, notebook_paths, page_entries, navigation_entries, display_title, registry_paths

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

"""Validate a non-notebook page path stored in the unified site registry."""
function normalize_page_entry(entry::AbstractString, filename::AbstractString; root::AbstractString=REPOSITORY_ROOT)
    relative_path = normpath(String(entry))
    isabspath(relative_path) && error("$(filename): absolute paths are not allowed: $(entry)")
    (relative_path == ".." || startswith(relative_path, ".." * Base.Filesystem.path_separator)) &&
        error("$(filename): paths must stay inside src: $(entry)")
    page_path = joinpath(root, "src", relative_path)
    isfile(page_path) || error("$(filename): page does not exist: src/$(relative_path)")
    relative_path
end

"""Use frontmatter's title when available, otherwise the first Markdown heading."""
function display_title(path::AbstractString, frontmatter=Dict())
    frontmatter_title = get(frontmatter, "title", nothing)
    frontmatter_title isa AbstractString && !isempty(strip(frontmatter_title)) && return strip(frontmatter_title)

    source = read(path, String)
    markdown_blocks = endswith(path, ".jl") ?
        (block.captures[1] for block in eachmatch(r"md\"\"\"(.*?)\"\"\""s, source)) :
        (source,)

    for block in markdown_blocks
        heading = match(r"(?m)^\s{0,3}#{1,6}\s+(.+?)\s*#*\s*$", block)
        heading === nothing && continue
        return strip(heading.captures[1])
    end

    splitext(basename(path))[1]
end

function load_registry(filename::AbstractString; root::AbstractString=REPOSITORY_ROOT)
    path = joinpath(root, filename)
    isfile(path) || error("Notebook registry not found: $(path)")
    YAML.load_file(path)
end

"""Read an explicitly ordered section list from the unified registry."""
function ordered_sections(registry, key::AbstractString, filename::AbstractString)
    sections = get(registry, key, nothing)
    sections isa AbstractVector || error("$(filename) must contain a '$(key):' list")

    result = NamedTuple{(:name, :data), Tuple{String, Any}}[]
    for section in sections
        section isa AbstractDict || error("$(filename): every $(key) entry must be an object")
        name = get(section, "name", nothing)
        name isa AbstractString && !isempty(strip(name)) || error("$(filename): every $(key) entry needs a non-empty 'name'")
        push!(result, (name=strip(name), data=section))
    end
    result
end

function registry_notebook_entries(filename::AbstractString; root::AbstractString=REPOSITORY_ROOT)
    registry = load_registry(filename; root)
    entries = NamedTuple{(:section, :kind, :path), Tuple{String, Symbol, String}}[]
    for section in ordered_sections(registry, "sections", filename)
        for kind in ("live", "static")
            listed = get(section.data, kind, nothing)
            listed === nothing && continue
            listed isa AbstractVector || error("$(filename): section $(repr(section.name)) field $(repr(kind)) must contain a list")
            for entry in listed
                entry isa AbstractString || error("$(filename): notebook entries must be strings")
                push!(entries, (section=section.name, kind=Symbol(kind), path=normalize_notebook_entry(entry, filename; root)))
            end
        end
    end

    paths = getproperty.(entries, :path)
    length(unique(paths)) == length(paths) || error("$(filename) contains a notebook more than once")
    entries
end

registry_paths(filename::AbstractString; root::AbstractString=REPOSITORY_ROOT) =
    getproperty.(registry_notebook_entries(filename; root), :path)

"""Return the unified registry entries, preserving the order within each YAML list."""
function notebook_entries(; root::AbstractString=REPOSITORY_ROOT)
    filename = "live-notebooks.yml"
    registry_notebook_entries(filename; root)
end

"""Return the YAML-defined non-notebook navigation pages, in display order."""
function page_entries(; root::AbstractString=REPOSITORY_ROOT)
    filename = "live-notebooks.yml"
    registry = load_registry(filename; root)

    entries = NamedTuple{(:section, :kind, :path), Tuple{String, Symbol, String}}[]
    for section in ordered_sections(registry, "pages", filename)
        listed = get(section.data, "entries", nothing)
        listed isa AbstractVector || error("$(filename): page section $(repr(section)) must contain a list")
        for entry in listed
            entry isa AbstractString || error("$(filename): page entries must be strings")
            push!(entries, (section=section.name, kind=:page, path=normalize_page_entry(entry, filename; root)))
        end
    end

    paths = getproperty.(entries, :path)
    length(unique(paths)) == length(paths) || error("$(filename) contains a page more than once")
    entries
end

"""Return every sidebar entry in the exact YAML-defined order."""
function navigation_entries(; root::AbstractString=REPOSITORY_ROOT)
    entries = vcat(page_entries(; root), notebook_entries(; root))
    paths = getproperty.(entries, :path)
    length(unique(paths)) == length(paths) || error("live-notebooks.yml contains a navigation entry more than once")
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
