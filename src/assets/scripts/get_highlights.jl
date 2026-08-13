let
    movies_dir = joinpath(@__DIR__, "..", "movies")
    gif_files = sort(filter(f -> lowercase(splitext(f)[2]) == ".gif", readdir(movies_dir)))

    highlight_lookup = Dict(
        basename(x["img"]) => x
        for x in get(metadata["homepage"], "highlights", [])
    )

    prettify(filename) = join(uppercasefirst.(split(replace(splitext(filename)[1], r"[_\-]+" => " "))), " ")

    slides = [
        let
            entry = get(highlight_lookup, gif, nothing)
            name = entry === nothing ? prettify(gif) : entry["name"]

            @htl("""
            <div class="carousel-slide gif-slide">
                <img src="$(root_url)/assets/movies/$(gif)" loading="lazy" alt="$(name)">
            </div>""")
        end
        for gif in gif_files
    ]

    isempty(slides) ? nothing : @htl("""
    <div class="carousel-section" id="explore">
      <div class="carousel" id="carousel-highlights">
        <button class="carousel-arrow carousel-prev" type="button" aria-label="Previous highlight">‹</button>
        <div class="carousel-track">
        $(slides)
        </div>
        <button class="carousel-arrow carousel-next" type="button" aria-label="Next highlight">›</button>
      </div>
    </div>
    """)
end
