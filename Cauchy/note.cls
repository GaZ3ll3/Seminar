%
% University of Texas at Austin,
% Department of Mathematics.
% 
% Template for Seminar note
% 
% Copyrigt@2014 Yimin Zhong
% 
% $id: note.cls release_version:1.0 By 2014/02/08 
%
% License(MIT)
%

\NeedsTeXFormat{LaTeX2e}[1999/12/01]
\ProvidesClass{note}[2014/02/08]
\def\version{1.0}
\LoadClass[cs4size,a4paper,fancyhdr,fntef]{article}

\RequirePackage{amsmath,amsthm,amsfonts,amssymb,bm,bbm}
% Reference:http://mirror.hmc.edu/ctan/macros/latex/contrib/bbm/bbm.pdf
\RequirePackage{graphicx}
\RequirePackage{enumerate}
\RequirePackage{multirow}
\RequirePackage{multicol}
\RequirePackage{subfig,epsfig}
\RequirePackage{xcolor}
\RequirePackage{float}
\RequirePackage{hyperref}
\RequirePackage{longtable}

\newtheorem{theorem}{THEOREM}[section]
\newtheorem{lemma}{LEMMA}[section]
\newtheorem{example}{EXAMPLE}[section]
\newtheorem{definition}{DEF}[section]
\newtheorem{axiom}{AXIOM}[section]
\newtheorem{property}{PROPERTY}[section]
\newtheorem{proposition}{PROPOSITION}[section]
\newtheorem{corollary}{COROLLARY}[section]
\newtheorem{remark}{REMARK}[section]
\newtheorem{condition}{CONDITION}[section]
\newtheorem{conclusion}{CONCLUSION}[section]
\newtheorem{assumption}{ASSUMPTION}[section]

\renewcommand\exp[1]{\mathrm{e}^{#1}}
\newcommand\dd{\mathrm{d}}
\newcommand\nn{\mathbf{n}}
\newcommand{\lp}{\triangle}
\newcommand{\Lp}{\nabla^2}
\newcommand{\bfgs}{\mathbf{BFGS}}
\newcommand*{\pfrac}[2]{\frac{\partial #1}{\partial #2}}
\newcommand{\wrt}{$\mathbf{w.r.t}\;$}

