### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ 6833c638-c7c8-4999-ae74-8ef5ceaf0569
using TikzPictures, LaTeXStrings

# ╔═╡ 2dca045e-df15-44ff-af48-b2eb240fdf92
tikz_default_options = raw"""
  background rectangle/.style={fill=white}, show background rectangle, 
  """

# ╔═╡ ec104bb8-4d21-11ee-186e-9d7c5bb7cb16
tikz_preamble = raw"""
  \usepackage{tikz}
\usepackage{booktabs} 
\usepackage{tikzducks}
  \usepackage{tikz-3dplot}
  \usetikzlibrary{fit, matrix, shapes.geometric}
  \tikzset{% use tikzset, not tikzstyle
      cell/.style={
          rectangle, rounded corners=5pt, draw,
      }
  }
  \tikzset{% use tikzset, not tikzstyle
      cellv/.style={
          rectangle, rounded corners=5pt, draw, rotate=90,
      }
  }
  \usepackage{xifthen}
  \usetikzlibrary{hobby}
  \usepackage{pgfplots}
  \usepackage{fontawesome}
  \usepackage{bm,amsfonts,amsmath}
  \usetikzlibrary{backgrounds,pgfplots.groupplots,snakes, math}
  \usepgfplotslibrary{patchplots}
  \pgfplotsset{try min ticks=2}
  \usepackage{pgfplotstable} 
  \usetikzlibrary{plotmarks,positioning,spy}
  \usetikzlibrary{shapes.geometric, arrows, fadings}
  \usepgfplotslibrary{groupplots, polar}
  \usepackage[space]{grffile}

  \usetikzlibrary{%
              decorations.pathreplacing,%
                  decorations.pathmorphing%
                  }
                  \usetikzlibrary{positioning,fit,backgrounds}


  \usetikzlibrary{backgrounds,
                hobby}
  \usetikzlibrary{shapes,arrows}
  \usetikzlibrary{decorations.markings}
  \usetikzlibrary{patterns}
  \usetikzlibrary{plotmarks}
  \usetikzlibrary{fit}
  \usetikzlibrary{intersections}
  \usepgfplotslibrary{fillbetween}

    \pgfplotsset{
                axis line style={black!10},
                    every axis label/.append style ={black!10},
                    every axis title/.append style ={black!10},
                        every tick label/.append style={black!10}  
                          }

  % need for pgfplots
  \newcommand{\axisz}{0cm}
  \newcommand{\axisx}{0cm}

  \usetikzlibrary{positioning}
  \usetikzlibrary{shapes.geometric}
  \usetikzlibrary{backgrounds}



\tikzset{
pics/wavelet/.style n args={3}{
  code={
\draw [fill=#2, thick, domain=-3:2, samples=200, smooth] plot coordinates{#1};
\node[fill=white] at (1.5,0) {#3};
   }
}}

% arguments (frequency, time shift, amplitude, color)

\tikzset{
pics/ricker/.style n args={4}{
  code={
\draw [fill=#4, thick, domain=-3:2, samples=200, smooth] plot (\x, {(1-(2*pi*(#1*(\x-#2))*(#1*(\x-#2))))* e^(-pi*(#1*(\x-#2))*(#1*(\x-#2))) * #3});
   }
}}


\tikzset{
pics/meyer/.style n args={4}{
  code={
\draw [fill=#4, thick, domain=-3:2, samples=200, smooth] plot (\x, {(sin(2*pi*#1*200*(\x-#2))-sin(pi*#1*200*(\x-#2))) / pi / (#1*200*(\x-#2)) * 20 * #3});
   }
}}

\tikzset{
pics/poisson/.style n args={4}{
  code={
\draw [fill=#4, thick, domain=-3:2, samples=200, smooth] plot (\x, {#3 * (1-(\x-#2)*(\x-#2)*#1*#1) / (1+(\x-#2)*(\x-#2)*#1*#1)/(1+(\x-#2)*(\x-#2)*#1*#1)});
   }
}}


\newcommand\waveletone[5]{
\begin{scope}[shift={#1}, rotate=0]
\draw [fill=#5, thick, domain=-3:2, samples=200, smooth] plot (\x, {(\x-#3)*(1-(2*pi*(#2*(\x-#3))*(#2*(\x-#3))))* e^(-pi*(#2*(\x-#3))*(#2*(\x-#3))) * #4});
\end{scope}
}

\newcommand\wavelettwo[5]{
\begin{scope}[shift={#1}, rotate=0]
\draw [fill=#5, thick, domain=-3:2, samples=200, smooth] plot (\x, {(1-(2*pi*(#2*(\x-#3))*(#2*(\x-#3))))* e^(-pi*(#2*(\x-#3))*(#2*(\x-#3))) * #4});
\end{scope}
}

\newcommand\waveletthree[5]{
\begin{scope}[shift={#1}, rotate=0]
\draw [fill=#5, thick, domain=-3:2, samples=200, smooth] plot (\x, {(\x*\x-#3)*(1-(2*pi*(#2*(\x-#3))*(#2*(\x-#3))))* e^(-pi*(#2*(\x-#3))*(#2*(\x-#3))) * #4});
\end{scope}
}


\newcommand\hx{6cm}
\newcommand\dx{0.5cm}
\newcommand\hz{3cm}
\newcommand\dz{0.5cm}


\tikzset{
pics/spike/.style={
  code={\draw[  draw=gray,
    thick,
    fill=gray,] (0,0) rectangle (0.2cm,0.2cm);
\draw[very thick,black] (.1cm,0) -- (.1cm,-0.2cm); }},
pics/receivers/.style n args={2}{
  code={
\foreach \x in {0,-0.8,-1.6,-2.4,-3.2,-4.0}{
\pic at (\x,0) {spike};
}
\draw[] (0,0) -- (-4,0);
\node[scale=2] at (-4.0,0.5cm) {receivers};
}
},
pics/truck/.style n args={6}{
  code={
\node[xshift=-2.6*\hx,yshift=#2*\hz] at (0,0) (s1) {};
\node[xshift=-\hx/4,yshift=#2*\hz] at (0,0) (s2) {};
\draw [very thick, black!90, <->,>=stealth] (s1)  -- (s2) ;
\filldraw[pattern=horizontal lines light gray	] (s1) rectangle ([yshift=-\hz/2]s2);
\draw[fill=gray!20] {(-#1*\hx-1.7cm,#2*\hz+0.4cm) -- (-#1*\hx+1.7cm,#2*\hz+0.4cm) -- (-#1*\hx+1.7cm,#2*\hz+2.2cm) -- (-#1*\hx+0.7cm,#2*\hz+2.2cm) -- (-#1*\hx+0.6cm,#2*\hz+1.6cm) -- (-#1*\hx-1.7cm,#2*\hz+1.6cm) -- (-#1*\hx-1.7cm,#2*\hz+0.4cm)};
\node at (-#1*\hx-1.7cm,#2*\hz+.8cm) (s9) {};
\draw[fill=gray!50] (-#1*\hx-1.3cm,#2*\hz+1cm) rectangle (-#1*\hx+0.3cm,#2*\hz+1.2cm);
\draw[fill=gray!5] (-#1*\hx+1cm,#2*\hz+1.2cm) rectangle (-#1*\hx+1.5cm,#2*\hz+1.7cm);
\draw[fill=gray!80] (-#1*\hx-1cm,#2*\hz+.4cm) circle (0.4cm) node(c1) {};
\draw[fill=gray!80] (-#1*\hx+1cm,#2*\hz+.4cm) circle (0.4cm);
\node[scale=2] at (-#1*\hx-.7cm,#2*\hz+1.9cm) {source};
\node at (-#1*\hx-3.9cm,#2*\hz) (s10) {};
\draw[thick] (s9.south) -- (s10.north);
\draw[snake=expanding waves,segment angle=90,#3,segment length=#4,thick] (-#1*\hx-1cm,#2*\hz) -- (-#1*\hx-1cm,#2*\hz-1.5cm);
\draw[snake=expanding waves,segment angle=90,#5,segment length=#6,thick] (-#1*\hx+1cm,#2*\hz) -- (-#1*\hx+1cm,#2*\hz-1.5cm);
\draw[snake=snake,#3,segment length=#4,thick] (-#1*\hx-2.7cm,#2*\hz-.5cm) -- (-2.55*\hx,#2*\hz-.5cm);
\draw[snake=snake,#5,segment length=#6,thick] (-#1*\hx-2.7cm,#2*\hz-.7cm) -- (-2.55*\hx,#2*\hz-.7cm);
\pic at (s10.north) {receivers={1}{2}};
   }
}}
\tikzset{
pics/crust/.style n args={6}{
  code={
 \coordinate (c0) at (-4, -1);
\coordinate (c1) at (-3,-1.1);
\coordinate (c2) at (-2,-1.2);
\coordinate (c3) at (-1,-1.3);
\coordinate (c4) at (0,-1.4);
\coordinate (c5) at (1,-1.6);
\coordinate (c6) at (2,-1.8);
	% Round rectangle
    \fill[gray!10] (-4,0) rectangle (3,-3);
    % Interface
    \draw[gray,line width=.5pt,draw,decorate,decoration={border,angle=-45,
                    amplitude=0.3cm,segment length=2mm}](-4,0)--(3,0);
    % Vertical dashed line
    \draw[gray!20](0,-3)--(0,0);
	\draw[black!20](-4,-3)--(3,-3);

    % Incidence

\foreach \x/\i/\y/\s/\c in {-2/1/-1.2/#1/#4,-0.5/2/-1.35/#2/#5,1/3/-1.6/#3/#6}{
\node[scale=1.5, star, star point ratio=0.5, star points=3, fill=black,anchor=center] at (\x,0) (rec\i) { };
\draw(rec\i)node[above]{receiver \#\i};
\draw[snake=snake,thick,segment length=\s,\c]
         (\x-0.2,-3)--(\x-0.1,\y)node[midway,right]{P};
\draw[snake=snake,thick,segment length=\s,\c, dashed](\x-0.1,\y)--(rec\i)node[midway,right]{S};
}
    % Media names
       \node[align=left]      at    (-3.4,-1.75) {subducting \\plate};
	
    % \draw[line width=.6pt] (0,0)
    %                       +(-135:.12cm) -- +(45:.12cm)
    %                       +(-45:.12cm) -- +(135:.12cm);
    % Interface pointer
    \draw[-latex,thick](2.25,0.5)node[right]{free surface}
         to[out=180,in=90] (2,0);

\draw[black, thick,pattern=horizontal lines,pattern color=gray] plot [smooth, tension=0.01] coordinates {(c0) (c1) (c2) (c3) (c4) (c5) (c6) (3,-1.9) (3,-2.3) (2,-2.1) (1,-1.9) (0,-1.8) (-1,-1.7) (-2,-1.6) (-3,-1.5) (-4,-1.5) };
   }}}
  """

# ╔═╡ 7035d9ef-d399-4339-a243-68d4a45bb425
plot(code; width="", preamble=raw"", options=raw"") = TikzPicture(code, options=tikz_default_options * options, preamble=tikz_preamble * preamble, width=width)

# ╔═╡ 5150711d-5c7d-4788-86d8-848cb072532b
tikz_default_options * raw"""scale=2"""

# ╔═╡ 78189b14-d836-4b54-bb26-1ed931f25521
plot(L"""
\draw[] (1,0) -- (0,1);
""", options=raw"""scale=2""")

# ╔═╡ c58e400e-44ed-40b1-bbed-0f73aecd1826
plot(L"""
\tikzset{
 pics/blob/.style={
   code={
   \draw[use Hobby shortcut, fill, closed] (0,0) +($(0:1+4*rnd)$)
       \foreach \a in {60,120,...,350} {  .. +($(\a: 1+4*rnd)$) };
   }
}}
  \pic at (0,0)  [fill=green!30, scale=1, rotate=360*rnd]{blob};
""")

# ╔═╡ f9e6fb55-852f-4d0a-aeb9-821b2471e500


# ╔═╡ 516564f8-2ccb-4d48-8435-0e8b486c30f6
md"# Ricker"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
TikzPictures = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"

[compat]
LaTeXStrings = "~1.3.1"
TikzPictures = "~3.5.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.2"
manifest_format = "2.0"
project_hash = "f5224af67e0622f73b77f865e6f1cd244c55def2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a2f1c8c668c8e3cb4cca4e57a8efdb09067bb3fd"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.0+2"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.0+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "7c82e6a6cd34e9d935e9aa4051b66c6ff3af59ba"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.2+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.HarfBuzz_ICU_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "HarfBuzz_jll", "ICU_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "6ccbc4fdf65c8197738c2d68cc55b74b19c97ac2"
uuid = "655565e8-fb53-5cb3-b0cd-aec1ca0647ea"
version = "2.8.1+0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.ICU_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "20b6765a3016e1fca0c9c93c80d50061b94218b7"
uuid = "a51ab1cf-af8e-5615-a023-bc2c838bba6b"
version = "69.1.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c84a835e1a09b289ffcd2271bf2a337bbdda6637"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "70c5da094887fd2cae843b8db33920bac4b6f07d"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "9fd170c4bbfd8b935fdc5f8b7aa33532c991a673"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.11+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fbb1f2bef882392312feb1ede3615ddc1e9b99ed"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.49.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg"]
git-tree-sha1 = "110897e7db2d6836be22c18bffd9422218ee6284"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.12.0+0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenJpeg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll", "LittleCMS_jll", "Pkg", "libpng_jll"]
git-tree-sha1 = "76374b6e7f632c130e78100b166e5a48464256f8"
uuid = "643b3616-a352-519d-856d-80112ee9badc"
version = "2.4.0+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a12e56c72edee3ce6b96667745e6cbbe5498f200"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.23+0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.Poppler_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "OpenJpeg_jll", "Pkg", "libpng_jll"]
git-tree-sha1 = "02148a0cb2532f22c0589ceb75c110e168fb3d1f"
uuid = "9c32591e-4766-534b-9725-b71a8799265b"
version = "21.9.0+0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TikzPictures]]
deps = ["LaTeXStrings", "Poppler_jll", "Requires", "tectonic_jll"]
git-tree-sha1 = "79e2d29b216ef24a0f4f905532b900dcf529aa06"
uuid = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"
version = "3.5.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "52ff2af32e591541550bd753c0da8b9bc92bb9d9"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.7+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e678132f07ddb5bfa46857f0d7620fb9be675d3b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d7015d2e18a5fd9a4f47de711837e980519781a4"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.tectonic_jll]]
deps = ["Artifacts", "Fontconfig_jll", "FreeType2_jll", "Graphite2_jll", "HarfBuzz_ICU_jll", "HarfBuzz_jll", "ICU_jll", "JLLWrappers", "Libdl", "OpenSSL_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "54867b00af20c70b52a1f9c00043864d8b926a21"
uuid = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"
version = "0.13.1+0"
"""

# ╔═╡ Cell order:
# ╠═6833c638-c7c8-4999-ae74-8ef5ceaf0569
# ╠═2dca045e-df15-44ff-af48-b2eb240fdf92
# ╠═ec104bb8-4d21-11ee-186e-9d7c5bb7cb16
# ╠═7035d9ef-d399-4339-a243-68d4a45bb425
# ╠═5150711d-5c7d-4788-86d8-848cb072532b
# ╠═78189b14-d836-4b54-bb26-1ed931f25521
# ╠═c58e400e-44ed-40b1-bbed-0f73aecd1826
# ╠═f9e6fb55-852f-4d0a-aeb9-821b2471e500
# ╟─516564f8-2ccb-4d48-8435-0e8b486c30f6
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
