# LaTeX Agent Guide

Read this guide before editing TeX, BibTeX, manuscripts, papers, or mathematical notes. These rules are about source style and notation discipline; broader build, citation, and manuscript policy will be added later.

## Source prose formatting

- Do not hard-wrap prose paragraphs in LaTeX source. Let the editor soft-wrap paragraphs visually.
- Preserve intentional blank lines between paragraphs, equations, and structural blocks.
- Do not reflow unrelated prose while making a local mathematical or formatting edit.

## Section markers

- Make sections easy to find in the LaTeX source by using prominent comment banners around `\section` and `\subsection` headings.
- Use this section template:

```tex
%%%%%%%%%%%%%%%%%%%%%%%%
%
\section{Title}
\label{sec:label}
%
%%%%%%%%%%%%%%%%%%%%%%%%
```

- Use this subsection template:

```tex
%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Title}
\label{subsec:label}
%%%%%%%%%%%%%%%%%%%%%%%%
```

- subsubsections should use the same format as subsections.

- chapters should use ths same format as sections.

- Preserve existing local label prefixes and naming conventions when adding labels.

## Equation formatting

- Do not add line breaks inside equations merely for source-code alignment.
- Keep an equation on one rendered line when it fits cleanly and compiles without overfull boxes.
- Use `aligned`, `align`, or manual equation line breaks only when needed for readability or page fit.
- When breaking an equation, keep related terms together and avoid splitting simple factors, bracketed products, short endpoint terms, or short operator applications across many lines.

## Notation discipline

- Minimize subscripts, superscripts, and indexes unless they distinguish genuinely different objects.
- Avoid redundant notation where an argument, surrounding context, or established symbol already carries the meaning.
- Prefer the document's established notation over inventing decorated variants.
- Do not introduce new notation just to make one local derivation feel more explicit.
