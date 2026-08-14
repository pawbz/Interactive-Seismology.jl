let
    function introductory_excerpt(notebook_path::AbstractString; limit::Int=180)
        source = read(notebook_path, String)
        for block in eachmatch(r"md\"\"\"(.*?)\"\"\""s, source)
            lines = filter(split(block.captures[1], '\n')) do line
                text = strip(line)
                !isempty(text) &&
                    !startswith(text, "#") &&
                    !startswith(text, "!!!") &&
                    !startswith(text, "Instructor:")
            end
            isempty(lines) && continue

            excerpt = replace(join(strip.(lines), " "), r"\s+" => " ")
            excerpt = replace(excerpt, r"\[([^\]]+)\]\([^\)]+\)" => s"\1")
            length(excerpt) > limit && return first(excerpt, limit - 1) * "…"
            return excerpt
        end
        nothing
    end

    slides = [
        let
            notebook_path = joinpath(@__DIR__, "..", "..", entry.path)

            name = entry.title
            desc = get(entry.frontmatter, "description", nothing)
            desc = desc isa AbstractString && !isempty(strip(desc)) ? desc : introductory_excerpt(notebook_path)
            image = get(entry.frontmatter, "image", nothing)
            href = root_url * "/" * splitext(entry.path)[1] * "/"

            @htl("""
            <a class="carousel-slide no-decoration" href="$(href)" title="$(desc)">
                <div class="slide-content">
                    <p class="slide-section">$(entry.section)</p>
                    <h3>$(name)</h3>$(desc === nothing ? nothing : @htl("<p>$(desc)</p>"))</div>$(image === nothing || isempty(image) ? nothing : @htl("""<div class="slide-preview"><img src="$(image)" loading="lazy" alt="$(name)"></div>"""))</a>""")
        end
        for entry in notebook_catalog
    ]

    isempty(slides) ? nothing : @htl("""
    <div class="carousel-section">
      <div class="carousel" id="carousel-notebooks">
        <button class="carousel-arrow carousel-prev" type="button" aria-label="Previous notebook">‹</button>
        <div class="carousel-track">
        $(slides)
        </div>
        <button class="carousel-arrow carousel-next" type="button" aria-label="Next notebook">›</button>
      </div>
    </div>
    """)
end
