# About

This page of the preferences tells which version of RetroShare you are
running and gives access to the credits and to the licence. It is the place to
look when you report a problem.

## Layout

- At the top, the buttons **About** and **Copy Info**.
- Below, a picture with the RetroShare logo and three lines:
  - **RetroShare version**: the version of the program. It is read as: the
    last release, then the number of changes made since that release, then
    the letter `g` followed by the identifier of the exact source revision,
    for instance `0.6.7-963-g49d11c2f`. A build made exactly from a release
    shows the release alone. A version that ends with `-dirty` was built from
    modified sources. "version not available" means the build did
    not record it.
  - **libretroshare version**: the version of the RetroShare engine, the
    library that does the network work behind the window.
  - **Qt version**: the version of the Qt toolkit the program was built with.

## Buttons

- **About**: opens the dialog described below.
- **Copy Info**: copies a technical description of your installation to the
  clipboard. It contains the RetroShare version, the libretroshare version,
  the name and version of your operating system, the Qt version, the list of
  the libraries used by RetroShare and by each loaded plugin with their
  versions, and the build options. Paste it into a bug report or a support
  request: it holds no personal data, no identifier of your node and no
  address.
- The picture: a click on it starts or stops an animation. A click with the
  Ctrl key held down swaps the picture for a small game; Ctrl+click again
  brings the picture back.

## Dialogs

### About dialog

The dialog opened by **About** has six tabs:

- **About**: the RetroShare version, a short presentation of the program and
  links to the web sites of the project, which open in your browser.
- **Authors**: the people who wrote RetroShare.
- **Thanks to**: the contributors and the projects RetroShare builds on.
- **Translation**: the translators, and where to help translate RetroShare.
- **License Agreement**: the text of the licence under which RetroShare is
  distributed.
- **Libraries**: the name and version of each library in use: first
  libretroshare and the libraries it uses, then one block per loaded plugin.
  The text can be selected and copied.

Close the dialog with the close button of its title bar.

## How to

### Give your version in a bug report

1. Click **Copy Info**.
2. Paste the text into the report, with the description of the problem.

## See also

- [Plugins preferences](settings-plugins.md)
- [Getting started](start.md)
