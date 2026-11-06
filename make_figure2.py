"""Builds the LaTeX/pgfplots source of Figure 2 (figure2.tex) directly from results.json."""
import json
r=json.load(open('results.json'))
G=[0.0,0.5,1.0,1.5,2.0]
def pts(n,key): return ''.join(f"({g},{r[f'moderate_{n}_{g}'][key]:.3f})" for g in G)
def t1(sc,n,key): return r[f'{sc}_{n}_0.0'][key]
bars=[]
# type I error bar chart: groups = moderate200, moderate400, strong200, strong400
grp=[('moderate',200),('moderate',400),('strong',200),('strong',400)]
def coords(key): return ' '.join(f"({i},{t1(sc,n,key):.3f})" for i,(sc,n) in enumerate(grp))
fig=r"""\begin{figure}[t]
\centering
\begin{minipage}[b]{0.42\linewidth}\centering
\begin{tikzpicture}
\begin{axis}[width=0.93\linewidth,height=4.1cm,font=\scriptsize,ybar,bar width=4.5pt,
  ymin=0,ymax=0.12,ylabel={type I error at $\alpha=0.05$},
  xtick={0,1,2,3},xticklabels={{moderate, 200},{moderate, 400},{strong, 200},{strong, 400}},
  x tick label style={align=center,text width=1.3cm},xlabel={state dependence, patients},
  enlarge x limits=0.18,grid=major,grid style={gray!20},
  ytick={0,0.05,0.10},yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=2},scaled y ticks=false,
  legend style={at={(0.5,1.02)},anchor=south,font=\tiny,draw=none,fill=none,legend columns=3,/tikz/every even column/.append style={column sep=3pt}}]
\addplot[fill=black!75,draw=black] coordinates {"""+coords('boot')+r"""};
\addplot[fill=gray!40,draw=black] coordinates {"""+coords('perm')+r"""};
\addplot[fill=white,draw=black] coordinates {"""+coords('lrt')+r"""};
\draw[red!70!black,dashed,thick] (axis cs:-0.5,0.05) -- (axis cs:3.5,0.05);
\legend{bootstrap,permutation,LRT}
\end{axis}
\end{tikzpicture}\\[-2pt]{\footnotesize\textbf{(a)} Calibration with no true configuration.}
\end{minipage}\hfill
\begin{minipage}[b]{0.52\linewidth}\centering
\begin{tikzpicture}
\begin{axis}[width=\linewidth,height=4.1cm,font=\scriptsize,
  xlabel={true configuration strength $\gamma$ (log odds ratio)},ylabel={rejection rate at $\alpha=0.05$},
  xmin=-0.1,xmax=2.1,ymin=0,ymax=1.02,grid=major,grid style={gray!20},xtick={0,0.5,1,1.5,2},
  legend style={at={(0.02,0.98)},anchor=north west,font=\tiny,fill=white,fill opacity=0.9,text opacity=1}]
\addplot[thick,black,mark=*,mark size=1.3pt] coordinates {"""+pts(200,'boot')+r"""};
\addplot[thick,red!75!black,mark=square*,mark size=1.3pt] coordinates {"""+pts(400,'boot')+r"""};
\addplot[thin,dashed,black,mark=o,mark size=1.2pt] coordinates {"""+pts(200,'lrt')+r"""};
\addplot[thin,dashed,red!75!black,mark=square,mark size=1.2pt] coordinates {"""+pts(400,'lrt')+r"""};
\addplot[black,dotted,thick,domain=-0.1:2.1] {0.05};
\legend{{bootstrap, $N=200$},{bootstrap, $N=400$},{in-sample LRT, $N=200$},{in-sample LRT, $N=400$},{$\alpha=0.05$}}
\end{axis}
\end{tikzpicture}\\[-2pt]{\footnotesize\textbf{(b)} Power of the bootstrap test, moderate dependence.}
\end{minipage}
\caption{\textbf{Planning simulation for the CIG test.} Each replicate drew $N$ patients with a binary tumor state (prevalence 0.4), a binary immune state whose prevalence depended on tumor state (0.6 versus 0.3 for moderate dependence, 0.85 versus 0.1 for strong dependence), and a binary stratum, then generated pCR from a logistic model with main effects of 0.7, a stratum effect of 0.3, and configuration strength $\gamma$. CIG was computed by five-fold held-out deviance. Each test used 199 null replicates for computational economy (the study uses 2,000), and each point averages 200 simulated cohorts. LRT, in-sample likelihood-ratio test, valid only under correct model specification. Code and results are available in the repository listed in Supporting Documentation.}
\end{figure}
"""
open('figure2.tex','w').write(fig)
print('fig2 written')
