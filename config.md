<!--
Add here global page variables to use throughout your website.
-->
+++
author = "Michael P. Friedlander"
mintoclevel = 2

github = "https://github.com/mpf/ubc-cpsc-406"

# Add here files or directories that should be ignored by Franklin, otherwise
# these files might be copied and, if markdown, processed by Franklin which
# you might not want. Indicate directories by ending the name with a `/`.
# Base files such as LICENSE.md and README.md are ignored by default.
ignore = ["node_modules/"]

# RSS (the website_{title, descr, url} must be defined to get RSS)
generate_rss = true
website_title = "UBC CPSC 406"
website_descr = "Computational Optimization"
website_url   = "https://friedlander.io/ubc-cpsc-406"
+++

\newcommand{\blurb}[1]{
    ~~~
    <span style="font-size:24px;font-weight:300;">!#1</span>
    ~~~
}
\newcommand{\refblank}[2]{
    ~~~
    <a href="!#2" target="_blank" rel="noopener noreferrer">#1</a>
    ~~~
}
\newcommand{\lineskip}{@@blank@@}
\newcommand{\skipline}{\lineskip}
\newcommand{\note}[1]{@@note @@title ⚠ Note@@ @@content #1 @@ @@}
\newcommand{\notation}[1]{@@note @@title Notation @@ @@content #1 @@ @@}
\newcommand{\warn}[1]{@@warning @@title ⚠ Warning!@@ @@content #1 @@ @@}
\newcommand{\theorem}[2]{ @@theorem **Theorem**: (_!#1_) #2 @@ }
\newcommand{\definition}[2]{ @@definition **Definition**: (_!#1_) #2 @@ }
\newcommand{\theorem}[2]{ @@theorem **Theorem**: (_!#1_) #2 @@ }
\newcommand{\algorithm}[2]{ @@theorem **Algorithm**: (_!#1_) #2 @@ }
<!--
Add here global latex commands to use throughout your pages.
-->
\newcommand{\T}{^T\!}
\newcommand{\R}{\mathbb R}
\newcommand{\ip}[1]{\langle #1 \rangle}
\newcommand{\minim}{\mathop{\hbox{\rm minimize}}}
\newcommand{\maxim}{\mathop{\hbox{\rm maximize}}}
\newcommand{\minimize}[1]{\displaystyle\minim_{#1}}
\newcommand{\maximize}[1]{\displaystyle\maxim_{#1}}
\newcommand{\st}{\hbox{\rm subject to}}
\newcommand{\range}{\hbox{range}}
\newcommand{\span}{\hbox{span}}
\newcommand{\Null}{\hbox{null}}
\newcommand{\norm}[1]{\|#1\|}

<!-- Matrices -->
\newcommand{\bmat}[1]{\begin{bmatrix}#1\end{bmatrix}}

\newcommand{\textt}[1]{\quad\text{#1}\quad}