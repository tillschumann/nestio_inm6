% ------------------------------------------------------------------------------------------------------------------------------------------------------
% YIC2015 class for LaTeX files                
% Based on the template of the YIC2012 and YIC2013
%
%
% The regular article class commands still exist and may be used.
% ------------------------------------------------------------------------------------------------------------------------------------------------------

% --------------------------------------------------
% Preliminary declaractions
% --------------------------------------------------

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{AbstractYIC2015}[2014/09/22]

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}} % to pass every unknown option to the article class
\ProcessOptions\relax % to process the options

\LoadClass[10pt,a4paper,twoside]{article}

% --------------------------------------------------
% Required packages 
% --------------------------------------------------
\usepackage{amsmath}
\usepackage{amsfonts}
\usepackage{amssymb}
\usepackage{times}

% --------------------------------------------------
% Page dimensions and definitions 
% --------------------------------------------------
\setlength\topmargin{0mm}
\setlength\headsep{13pt}
\setlength\oddsidemargin{0.0mm}
\setlength\evensidemargin{0.0mm}
\setlength\textwidth{162mm}
\setlength\textheight{235mm}

% --------------------------------------------------
% Frontmatter (title page) definitions 
% --------------------------------------------------
\renewcommand{\maketitle}{%
	\thispagestyle{empty}
		% --------------------------------------------------
		% Top Minipage environment definitions 
		% --------------------------------------------------
		\begin{minipage}[t]{\textwidth}
		\begin{flushright}
		\scriptsize
		% ---------------------------- PREPRINT HEADER -----
		\textsf{\textbf{YIC GACM 2015}} \\
		\textsf{3\textsuperscript{rd} ECCOMAS Young Investigators Conference} \\%
		\textsf{6\textsuperscript{th} GACM Colloquium} \\
		\textsf{July 20--23, 2015, Aachen, Germany} \\
		\end{flushright}
		\end{minipage}
		
		% --------------------------------------------------
		% Minipage environment definitions 
		% --------------------------------------------------
		\begin{minipage}[t]{\textwidth}
		{\raggedright
		% ---------------------------- TITLE -----
			\vspace{40pt}
				\fontsize{14}{16}\selectfont
				\bfseries
				\@title
				\par
		% ---------------------------- AUTHORS -----
			\vspace{8pt}%
				\fontsize{12}{14}\selectfont
				\normalfont
				\@author
				\par
		% ---------------------------- AFFILIATIONS -----
			\vspace{16pt}%
				\fontsize{8}{10}\selectfont
				\normalfont\@address
				\par
		% ---------------------------- CORRESPONDING AUTHOR -----
			         \fontsize{8}{10}\selectfont
			         \normalfont$^*$\@corauth
			         %\par  
				\vspace{24pt}%
		}	
			\par
		% ---------------------------- KEYWORDS -----
			\vspace{12pt}%
			\noindent\textbf{Keywords:\space}{\@keywords}
			\par
		% ---------------------------- ABSTRACT -----
			\vspace{16pt}
			\@abstract
			\par
			\vspace{4pt}
		\end{minipage}
}
% End of \@maketitle

% --------------------------------------------------
% New commands and definitions 
% --------------------------------------------------
\def\author{\@ifnextchar[{\@getshortauthor}{\@getauthor}}
\def\@getauthor#1{\gdef\@author{#1}\gdef\@shortauthor{#1}}
\def\@getshortauthor[#1]#2{\gdef\@shortauthor{#1}\gdef\@author{#2}}

\newcommand{\address}[1]{\gdef\@address{#1}}
\newcommand{\dates}[1]{\gdef\@dates{#1}}
\newcommand{\corauth}[1]{\gdef\@corauth{#1}}
\newcommand{\keywords}[1]{\gdef\@keywords{#1}}
\renewcommand{\abstract}[1]{\gdef\@abstract{#1}}
\newcommand{\heading}[1]{\gdef\@heading{#1}}

\newcommand{\authref}[1]{\,\ensuremath{^{\textrm{#1}}}}
\newcommand{\corref}[0]{\ensuremath{^{,*}}}
\newcommand{\authaddr}[2]{\authref{#1} #2\\ \vspace{6pt}}
			
% --------------------------------------------------
% The bibliography environment
% --------------------------------------------------
\renewenvironment{thebibliography}[1]
{
	\subsubsection*{\refname}%
		\list{\@biblabel{\@arabic\c@enumiv}}%
		{\settowidth\labelwidth{\@biblabel{#1}}%
			\leftmargin\labelwidth
				\advance\leftmargin\labelsep
				\@openbib@code
				\usecounter{enumiv}%
				\let\p@enumiv\@empty
				\renewcommand\theenumiv{\@arabic\c@enumiv}}%
				\sloppy
				\clubpenalty4000
				\@clubpenalty \clubpenalty
				\widowpenalty4000%
				\sfcode`\.\@m
				\topsep=0.0ex
				\parsep=0.0ex
				\parskip=0.0ex
				\itemsep=0.0ex
				\small

}
{\def\@noitemerr
	{\@latex@warning{Empty `thebibliography' environment}}%
	\endlist}
	
	% --------------------------------------------------
	% Space and intentation definitions
	% --------------------------------------------------
	\setlength{\textfloatsep}{12pt}
	\setlength{\abovedisplayskip}{6pt}
	\setlength{\belowdisplayskip}{6pt}
	\setlength{\parindent}{0mm}
	\setlength{\abovecaptionskip}{6pt}
	\setlength{\belowcaptionskip}{6pt}
	\setlength{\itemsep}{0pt}
	%\addtolength{\parskip}{2mm} %gilac
	
	% --------------------------------------------------
	% Final commands
	% --------------------------------------------------
	\AtBeginDocument{%
		\pagestyle{myheadings}
			\maketitle
			\let\maketitle\relax
			\flushbottom
	}

\endinput