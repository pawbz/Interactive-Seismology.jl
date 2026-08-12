if isempty(get(metadata["homepage"], "highlights", []))
    nothing
else
    highlights = [
      @htl("""
      <section>
      <div class="content">
          <h2>$(x["name"])</h2>
          <p>$(x["text"])</p>
      </div>
      <div class="preview">
          <img src="$(x["img"])">
      </div>
      </section>
      """) for x in metadata["homepage"]["highlights"]
    ]

    @htl("""
    <div id="explore" class="subjectscontainer wide">
      <div class="section-heading">
        <p class="section-kicker">Explore by experiment</p>
        <h2>Follow the signal</h2>
        <p>Begin anywhere—each notebook is a self-contained experiment with live controls and visual feedback.</p>
      </div>
      <div class="contain">
      $(highlights)      
      </div>
    </div>
    """)
end
