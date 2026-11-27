# Files

The Files page, titled **File Transfers**, is where you share files, find the
files of other people, and follow your downloads and uploads.

RetroShare transfers files in two ways: directly from your friends, and
anonymously from people further away in the network, through tunnels that
cross several nodes. A file can come from several sources at once, and you
become a source for the parts you already have while you are still
downloading.

## Layout

The page is a set of tabs. The tab that was open when you closed RetroShare is
open again the next time.

- **Downloads**: the list of your downloads at the top and, below a movable
  separator, the **Uploads** list of what other people are fetching from you.
  The number of entries is shown in each title. See
  [Downloads and uploads](files-transfers.md).
- **Search**: searches the files of your friends and of the rest of the
  network by keyword. See [Search](files-search.md).
- **Friends files**: the shared directories that your friends let you browse.
  See [Shared files](files-shared.md).
- **My files**: the directories you share, and the tools to configure them.
  See [Shared files](files-shared.md).

A plugin can add a tab of its own after these four.

## Buttons

- **?** (help button, top right): shows a short introduction to the page.
  Click the introduction to close it.

The page has no other button of its own: every action is in the tabs and in
their context menus.

## How file sharing works

### Shared directories

You share files by declaring directories in the share manager (**My files**
tab, **Configure shared directories** button). RetroShare reads every file of
these directories and computes its hash. Nothing is copied or uploaded at that
point: the files stay where they are, and they are sent only to the nodes that
ask for them.

Each shared directory has three independent permissions:

- **Browsable**: your friends, or only the groups of friends you choose, see
  the directory and its content in their **Friends files** tab, and find it
  with the **Browsable files** search.
- **Anonymous search**: the files of the directory answer the **Distant**
  searches made by anybody in the network. The person who searches does not
  learn who you are.
- **Anonymous download**: the files can be downloaded through anonymous
  tunnels by anybody who knows their hash.

A directory that allows anonymous download but not anonymous search holds
files that only the people who already have a link to them can fetch.

### Hashes

A file is identified by its hash, a 40 character fingerprint computed from its
content, not by its name. Two files with the same content are the same file
for RetroShare, whatever they are called, and each of their holders is one
more source. A new or modified file cannot be searched, linked or downloaded
until RetroShare has finished hashing it.

### Friends and tunnels

When you start a download, RetroShare asks your friends that have the file,
and also looks for it in the rest of the network. A source found beyond your
friends is reached through a tunnel: the data goes from friend to friend, and
none of the nodes on the way knows both ends. In the transfer lists such a
source is named **Tunnel via** followed by the friend through which the tunnel
goes.

### File links

A RetroShare file link holds the name, the size and the hash of a file. A link
can also describe a whole directory with all its files. You copy a link with
the **Copy RetroShare Link** entries of the context menus and you paste it in
a chat, a mail, a forum or a channel post. Whoever clicks it, or pastes it in
the **Downloads** list, can then download the file. See
[RetroShare links](links.md).

### Collections

A collection is a small file with the extension `.rscollection` that lists
files, organised in directories, with their sizes and hashes. It does not
contain the files themselves. You send it or share it like any other file;
whoever opens it chooses which of the listed files to download.

The **Collection** submenu is in the context menus of the **Downloads** list,
of the search results and of the shared files:

- **Create Collection...**: asks for the name of the collection file, writes
  the selected files into it and opens the collection editor.
- **Modify Collection...**: opens the selected `.rscollection` file in the
  collection editor. Available when the selection is one collection file that
  you have on your disk.
- **View Collection...**: opens the selected `.rscollection` file in the
  download window, where you choose what to download. Same condition.
- **Download from collection file...**: when the selection is one collection
  file that you have finished downloading, opens it in the download window.
  Otherwise asks you to pick a `.rscollection` file on your disk.

## Dialogs

### Collection

The window is titled **Collection** followed by the name of the file. It has
two forms.

**Collection Editor**, to build or change a collection:

- **File name :** is the collection file. The **...** button next to it
  chooses another file. When that file already exists you are asked to
  **Overwrite** it, **Merge** it or **Cancel**.
- **Total size :** and **Selected files :** sum up the content.
- The left pane is the tree of the files of your computer. The right pane is
  the content of the collection, with the columns **File**, **Size**,
  **Hash** and **Count** (the number of files in a directory).
- Between the two panes, four buttons:
  - add the item selected on the left, without the content of a directory;
  - add the item selected on the left with everything it contains;
  - remove the item selected on the right;
  - create a directory inside the collection. You are asked for its name.
    Select at most one directory of the collection first: the new one is
    created inside it.
- Added files must be hashed. Until it is done they show
  **[File is being hashed]** and **Save** is disabled.
- **Save**: writes the collection file and closes the window.
- **Cancel**: closes the window without saving.

**Download files**, to fetch the content of a collection:

- The tree shows the content of the collection with a check box on every
  line. Checking a directory checks everything inside it. **Selected files :**
  counts the checked files.
- **Destination:** is the directory where the files will be saved, by default
  your incoming directory. The **...** button, or a right click on the field,
  opens a menu with your shared directories that are writable and
  **Specify...**, which lets you pick any directory. The directories of the
  collection are created below the destination.
- **Download!**: starts the download of the checked files and closes the
  window. If a file name contains characters that your system does not accept,
  you choose between **Correct filename** (the characters are replaced by
  `_`), **Correct all**, **Skip this file** and **Skip all**.
- **Cancel**: closes the window without downloading.

## How to

### Share a directory

1. Open the **My files** tab and click **Configure shared directories**.
2. Click **Add new** and choose the directory.
3. Set its permissions with the three buttons of the **Access** column.
4. Click **Apply and close**.

### Send a file to someone

1. Put the file in a shared directory and wait until it is hashed.
2. In **My files**, right click the file and choose
   **Copy retroshare Links to Clipboard**.
3. Paste the link in a chat or a mail.

### Download from a link

1. Copy the link.
2. Open the **Downloads** tab, right click in the list and choose
   **Paste RetroShare Link**.
3. Check the files and the destination, then click **Download!**.

## Related settings

- [Files](settings-files.md): shared directories, incoming and partials
  directories, automatic check of the shared directories, number of
  simultaneous downloads, default chunk strategy, direct downloads.

## See also

- [Downloads and uploads](files-transfers.md)
- [Search](files-search.md)
- [Shared files](files-shared.md)
- [RetroShare links](links.md)
- [Statistics](statistics.md)
- [Glossary](glossary.md)
