# The RetroShare manual

This directory holds the manual that the GUI shows in its help browser
(F1, the "Help" entry of the tray menu, the "Read more in the manual" link of
the help bubbles).

- `en/contents.xml` is the table of contents. Each `Topic` has an `id`, the
  `name` shown in the tree and the Markdown file of the page (`html`
  attribute). Topic ids are joined with dots: `settings.network`.
- `en/*.md` are the pages. Every page must be listed in `content.qrc`.
- A translation goes in a directory named after the language code (`fr/`,
  `de/`...), with the same file names. A page that is not translated is shown
  in English.
- `tools/help-coverage.py` lists the pages and menu entries of the GUI that the
  manual does not mention. Run it after a change of the GUI.

A page of the GUI is tied to its topic in `MainWindow::initStackedPage()` (main
pages) and `SettingsPage::initStackedWidget()` (preferences).

## Writing a page

The reader is a user of RetroShare, not a developer: no class names, no file
names, no source code.

### Markdown subset

The pages are rendered by Qt and must stay plain CommonMark:

- headings, paragraphs, bullet and numbered lists (two levels at most), bold,
  inline code, links;
- no tables, no raw HTML, no images, no footnotes;
- a link to another page is the file name: `[Network](network.md)`.

### Structure

One page per screen. The sections come in this order; a section that has
nothing to say is left out.

    # Name of the screen, as the GUI shows it

    One or two paragraphs: what the screen is for.

    ## Layout
    The areas of the screen, from top to bottom and left to right.

    ## Buttons
    One bullet per button or control that is not self-explanatory.

    ## Context menus
    ### Where the right click happens
    One bullet per entry, in the order of the menu.

    ## Dialogs
    ### Title of the dialog
    The windows opened from this screen, with their fields.

    ## How to
    Short numbered recipes for the common tasks.

    ## Related settings
    Links to the pages of the preferences that change this screen.

    ## See also

### Rules

- Every button, option and menu entry of the screen is described.
- A label is written in bold, exactly as the English GUI shows it, without the
  `&` of the keyboard accelerators: `- **Copy RetroShare link**: copies...`.
  An entry of a submenu is written `**Parent > Child**`.
- An entry that is shown only in some cases says so ("only for your own
  identities", "not on a hidden node").
- The text says what the entry does, as the code does it. Nothing is guessed.
- Features that are compiled out are not documented.
- Second person, present tense, short sentences. The product name is spelled
  "RetroShare".
