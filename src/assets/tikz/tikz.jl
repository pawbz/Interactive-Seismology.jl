### A Pluto.jl notebook ###
# v0.20.19

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
