Accessibility
=============

This page states how accessible this documentation site is, what was checked
to say so, and what is known not to meet the target.

What this site targets
----------------------

The target is the Web Content Accessibility Guidelines (WCAG) 2.1, level AA.
The review is dated **2026-10-08**, and covered the site as built with
Sphinx 9.1.0, the alabaster 1.0.0 theme and nbsphinx 0.9.8. The measurements
are from the built HTML and CSS, not from a browser or a screen reader: no
assistive technology was run, so this is a conformance review of the markup
and the colours, not a user test. The page language is declared as English
(``<html lang="en">``, set in ``conf.py``).

What was checked
----------------

**Contrast.** The ratio is the WCAG relative-luminance ratio,
``(L1 + 0.05) / (L2 + 0.05)``. Normal text needs 4.5, large text and
non-text parts of the interface need 3. Every colour that renders text in the
theme is on a white page unless stated.

.. list-table::
   :header-rows: 1
   :widths: 34 16 12 10 10 18

   * - Element
     - Colours
     - Ratio
     - Needs
     - Result
     - After override
   * - Body text
     - ``#3E4349``
     - 9.98
     - 4.5
     - pass
     -
   * - Links
     - ``#004B6B``
     - 9.46
     - 4.5
     - pass
     -
   * - Links on hover
     - ``#6D4100``
     - 8.73
     - 4.5
     - pass
     -
   * - Sidebar text, links, lists
     - ``#555``, ``#444``, ``#000``
     - 7.46, 9.74, 21.0
     - 4.5
     - pass
     -
   * - Inline code
     - ``#222`` on ``#ecf0f3``
     - 13.88
     - 4.5
     - pass
     -
   * - Narrow-screen sidebar text and links
     - ``#FFF``, ``#AAA`` on ``#333``
     - 12.63, 5.44
     - 4.5
     - pass
     -
   * - Footer text and links
     - ``#888``
     - 3.54
     - 4.5
     - **fail**
     - ``#595959``, 7.00
   * - Heading permalink (the pilcrow)
     - ``#DDD``
     - 1.36
     - 3
     - **fail**
     - not changed (see below)
   * - Sidebar rule
     - ``#AAA``
     - 2.32
     - 3
     - decorative
     -

Code blocks have no background of their own in this theme, so the syntax
colours were measured on white. The default style's token colours that fail
were darkened to the nearest shade that passes:

.. list-table::
   :header-rows: 1
   :widths: 40 16 12 16 16

   * - Token (pygments class)
     - Colour
     - Ratio
     - Result
     - After override
   * - Comments (``.c``, ``.sd``)
     - ``#8F5902``
     - 5.84
     - pass
     -
   * - Keywords and builtins (``.k``, ``.nb``, ``.ow``)
     - ``#004461``
     - 10.50
     - pass
     -
   * - Numbers (``.m``)
     - ``#900``
     - 8.92
     - pass
     -
   * - Operators (``.o``)
     - ``#582800``
     - 12.22
     - pass
     -
   * - Strings (``.s``, ``.se``, ``.si``)
     - ``#4E9A06``
     - 3.53
     - **fail**
     - ``#3A7A04``, 5.30
   * - Console output, decorators (``.go``, ``.nd``)
     - ``#888``
     - 3.54
     - **fail**
     - ``#595959``, 7.00
   * - Attributes (``.na``)
     - ``#C4A000``
     - 2.51
     - **fail**
     - ``#8A6D00``, 4.92
   * - Labels, entities (``.nl``, ``.ni``)
     - ``#F57900``, ``#CE5C00``
     - 2.76, 4.07
     - **fail**
     - ``#A04800``, 6.14

The overrides are in ``docs/source/_static/custom.css``, one comment each
naming the measured ratio. Nothing that passed was restyled.

**Keyboard.** The built ``index.html`` and ``reading_the_diagnostics.html``
were read for their focusable elements in document order.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Check
     - Finding
   * - Focus order
     - The page content (96 links on the home page, 29 on
       *Reading the diagnostics*), then the sidebar (19 and 25 links and the
       search form), then the footer (3 links). No ``tabindex`` is set.
   * - Search
     - A native ``<form>`` with a text input (labelled through
       ``aria-labelledby``) and a submit button, reachable by Tab and
       submitted with Enter. No scripting is needed to use it.
   * - Visible focus
     - The theme's stylesheets set no ``outline: none``, so the browser's own
       focus indicator applies; no focus style was added.
   * - Skip link
     - There is none. The content comes first in the document, so none is
       needed to reach it, but see the limits below.

**Figures.** Every ``figure``, ``image`` and ``plot`` directive in the
reStructuredText pages carries ``:alt:``: the nine figures of
:doc:`reading_the_diagnostics` all do, and each describes what the figure
shows. The notebook pages are the exception, below. The palettes of the
library's plotting functions were checked for colour-vision deficiency, and
the record, with the one docs figure whose red and green lines merge, is on
the results API page (:doc:`ampere.results`).

What fails or is limited
------------------------

* **Notebook output images have no useful alt text.** nbsphinx writes the
  image file name as the alt text of each figure a notebook cell draws,
  because a matplotlib inline figure carries no alt metadata and the
  notebooks are executed when the site is built, so none can be stored in
  them. Two remedies exist: pass alt text in every plotting cell's output
  metadata, or add a Sphinx transform that supplies it. Neither is done yet;
  the convention check skips notebook pages because it could not pass them.
* **The search box is late in the tab order.** The theme puts the page's
  content before the sidebar, so reaching search by keyboard takes every
  link on the page first. This is a property of alabaster; the site does not
  override its templates.
* **Heading permalinks cannot be reached by keyboard** (they are hidden
  until the pointer is over the heading, and their colour is below 3:1).
  They duplicate the table of contents, so nothing is lost, and they were
  left as the theme has them.
* **Library figures in a few cases rely on hue alone** (traces, coverage
  curves, several anomaly scores on one axes). They are listed, with the
  remedy, on the results API page.

The figure convention
---------------------

Every ``figure``, ``image`` or ``plot`` directive in a ``.rst`` page must
carry an ``:alt:`` option that says what the figure *shows* (the finding a
reader should take from it), not what it is. A small extension in
``docs/source/conf.py`` warns about any that does not, and because the
documentation is built with warnings as errors, the build then fails.

Reporting a problem
-------------------

If something here is hard to use, please open a
`Documentation gap <https://github.com/ICSM/ampere/issues/new?template=docs.yml>`_
issue, naming the page and what got in the way. Questions go to
`Discussions <https://github.com/ICSM/ampere/discussions>`_.
