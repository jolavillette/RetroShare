# Shared files

Two tabs of the [Files](files.md) page show shared files. **My files** lists
the directories you share and is where you configure them. **Friends files**
lists the directories that your friends let you browse, and lets you download
from them. Both tabs look and work the same way; the differences are noted
below.

## Layout

- A bar at the top holds the buttons and the filter field described under
  Buttons.
- The rest of the tab is the list of files.

In **My files** the list has two top entries: **My files**, which holds your
shared directories under their visible names, and
**Temporary shared files**, which holds the single files that RetroShare
shares for a limited time, for instance the files you attached to a mail, a
chat message or a channel post.

In **Friends files** the list has one entry per friend that shares something
you are allowed to see. A friend shows up here only if they made at least one
directory browsable for you. Hovering a friend whose list has not arrived yet
shows **This node hasn't sent any directory information yet.** Files that you
already have, or that you are downloading, are shown in a different colour.

While RetroShare refreshes a list, the header of the first column ends with
**[updating...]**. A file that is not hashed yet has a special icon; it
cannot be linked or downloaded until hashing is done.

### Columns of My files

In tree view:

- **My Directories**: the name of the directory or the file.
- **# Files**: the number of files in the directory, sub-directories included.
- **Size**: the size of the file, or the total size of the directory.
- **Age**: the time since the file was last modified; for a directory, the
  age of its most recent file.
- **Access**: for a directory, three icons that recall its permissions:
  anonymous search, anonymous download, browsable by friends.
- **Visibility**: for a browsable directory, who can browse it:
  **[All friend nodes]**, or **Only** followed by the groups of friends you
  chose.
- **Uploaded**: the total amount of data that has been uploaded from this
  file, or from all the files of the directory, as far as RetroShare kept the
  statistics.

In flat view **# Files** and **Access** are hidden, and the column after
**Age** is **Directory**, the directory that holds the file.

### Columns of Friends files

In tree view: **Friends Directories**, **# Files**, **Size** and **Age**.

In flat view: **Friends Directories** (the file name), **Size**, **Age**,
**Friend** (the friend that shares the file) and **Directory**.

## Buttons

- **Configure shared directories** (My files only): opens the share manager
  described below.
- Filter field (magnifying glass icon): searches the list by name. Type at
  least 3 characters, then press Enter or click the button that appears at
  the right of the field (tooltip **Start Search**). With several words, only
  the names that contain all of them are kept. The list is unfolded to show
  the matches, and the tooltip of the field gives their number. The other
  button (tooltip **Reset**) empties the field and shows everything again.
  At most 3000 matches are handled.
- View list:
  - **Tree view**: directories and sub-directories, as on the disk.
  - **Flat view**: all the files in a single list without directories, handy
    to sort by size or age. The flat list has a maximum length; when it is
    reached, a tooltip on the list tells you that no more files are added.
- **Popular files** (My files only; tooltip **Show only files and folders
  that have been successfully uploaded to other peers**): keeps only the
  files that have been uploaded at least once, and the directories that
  contain some. In tree view the list is unfolded.
- Age list: **All**, **One day old**, **One Week old**, **One month old**.
  With one of the three periods, the list is sorted by **Age**, the folder
  icons change colour, and in the columns where it is displayed the files
  modified within the period are marked **NEW**. **All** goes back to the
  sort by name.
- **Check files** (My files only): forces a complete check of all your shared
  directories. The automatic check only looks for new and removed files;
  this button also hashes again the files that have changed. The button reads
  **Checking...** while a check is running.
- **Download** (Friends files only; tooltip **Download selected**): downloads
  the selected files and directories at once into your incoming directory.
  A directory is recreated there with its sub-directories.

The view, the sort order and the width of the columns are remembered for each
tab.

## Context menus

### My files: a shared directory

- **Open Folder**: opens the directory in the file manager of your system.
- **Copy retroshare Links to Clipboard**: copies a link that describes the
  whole directory, with all its files. Files that are not hashed yet are left
  out, and a message tells you so.
- **Retroshare Collection > Create Collection...**: creates a collection file
  with the content of the selected directories. See [Files](files.md).
- **Share on channel...**: a submenu with the channels in which you can
  publish. Choosing one opens a new post for that channel with the selected
  files attached. Only files are attached: select the files, not the
  directory that holds them. See [Channels](channels.md).
- **Share on forum...**: a submenu with the forums you are subscribed to.
  Choosing one opens a new message for that forum that contains the links of
  the selection. See [Forums](forums.md).

### My files: a file of a shared directory

- **Open File**: opens the file with the default application of your system.
- **Copy retroshare Links to Clipboard**: copies one file link per selected
  file.
- **Send retroshare Links**: opens a new mail that contains the links of the
  selection. See [Mail](mail.md).
- **Retroshare Collection > Create Collection...**: creates a collection file
  that lists the selected files.
- **Retroshare Collection > Modify Collection...**: shown only for a
  `.rscollection` file. Opens it in the collection editor.
- **Recommend in a message to...**: opens a new mail with the selected files
  listed as recommended files. The person who receives it can download them
  from the mail.
- **Share on channel...** and **Share on forum...**: as for a directory.

### My files: a temporary shared file

- **Open File**: opens the file.
- **Open Folder**: opens the directory that holds the file.
- **Copy retroshare Links to Clipboard**: copies the file link.
- **Send retroshare Links**: opens a new mail that contains the link.
- **Stop sharing this file**: removes the file from the temporary shared
  files. The file itself is not deleted.

The two top entries, **My files** and **Temporary shared files**, have no
menu.

### Friends files: a file or a directory

- **Download**: shown when a single file is selected. Starts the download at
  once. It appears in the [Downloads](files-transfers.md) tab.
- **Download...**: shown instead of **Download** for a directory, or when
  several lines are selected. Opens the **Download files** window, where you
  check what you want and choose the destination. See [Files](files.md).
- **Copy retroshare Links to Clipboard**: copies the links of the selection: a
  file link per file, a directory link per directory.
- **Send retroshare Links**: opens a new mail that contains these links.
- **Recommend in a message to...**: opens a new mail with the selected files
  listed as recommended files.
- **Collection > Create Collection...**: creates a collection file that lists
  the selection.
- **Collection > Download from collection file...**: asks for a collection
  file on your disk and opens it to download its content.
- **Show Empty** (tree view only, checkable): also lists the friends that
  share nothing you can see.

On a friend's name, or on the empty part of the list, the menu holds only
**Show Empty**, in tree view.

## Dialogs

### RetroShare Share Manager

Opened by **Configure shared directories**, and from the
[Files preferences](settings-files.md). It lists your shared directories, one
per line. What you change here is applied only when you click
**Apply and close**.

The columns:

- **Shared directory**: the directory on your disk. Double click it to choose
  another directory. It is greyed when the directory no longer exists.
- **Visible name**: the name under which your friends see the directory,
  instead of its real path. Double click it to type a name. A new directory
  shows **[Unset] (Double click to change)** until you name it. Two
  directories cannot have the same visible name.
- **Access**: three buttons that you switch on and off. Their tooltips tell
  the current state.
  - Anonymous search: when on, **Directory can be searched anonymously**: the
    files answer the distant searches of anybody in the network.
  - Anonymous download: when on, the files can be downloaded through
    anonymous tunnels by anybody who has their hash. If anonymous search is
    off for the same directory, the tunnels are also end-to-end encrypted.
    When off, **Files cannot be downloaded anonymously**.
  - Browsable: when on, **Directory content is visible to friend nodes (see
    list at right)**: the friends listed in the next column see the directory
    in their **Friends files** tab.
- **Visibility**: the friends that can browse the directory:
  **[All friend nodes]**, or the groups you chose. Double click it to open a
  window where you check groups of friends; with no group checked, all your
  friends can browse. The text is greyed while the directory is not
  browsable, because the choice has no effect then.

A directory you add starts with anonymous search and anonymous download on,
and browsable off: its files can be found and fetched anonymously, but your
friends do not see the directory in their lists.

The buttons:

- **Add new** (tooltip **Add a Share Directory**): asks for a directory and
  adds it to the list. A directory that is already in the list is not added
  twice.
- **Apply and close**: saves the list and closes the window. RetroShare then
  starts hashing the new files.
- **Cancel**: closes the window and drops your changes.

Right click a cell of the list. The first entry depends on the column:

- **Choose directory to share...** (on **Shared directory**): same as the
  double click.
- **Choose visible name...** (on **Visible name**): edits the name.
- **Change group visibility...** (on **Visibility**): opens the choice of
  groups.
- **Remove**: after a confirmation, takes the directory out of the list. The
  directory and its files stay on your disk; they are just no longer shared.

There is no menu on the **Access** column.

## How to

### Share a directory with your friends only

1. In **My files**, click **Configure shared directories**, then **Add new**,
   and choose the directory.
2. In the **Access** column, switch anonymous search and anonymous download
   off and browsable on.
3. To restrict it to some friends, double click the **Visibility** cell and
   check the groups.
4. Click **Apply and close**.

### Download a directory from a friend

1. Open **Friends files** and unfold the friend.
2. Right click the directory and choose **Download...**.
3. Uncheck what you do not want, choose the destination and click
   **Download!**.

### See what your friends added recently

1. Open **Friends files** and choose **Flat view**.
2. Choose **One Week old** in the age list: the most recent files come first.

### Make RetroShare notice a modified file

1. Open **My files**.
2. Click **Check files** and wait until the button no longer reads
   **Checking...**.

## Related settings

- [Files](settings-files.md): **Incoming Directory**,
  **Automatically share incoming directory (Recommended)**, the automatic
  check of the shared directories and its period, symbolic links, ignored
  file names, how long the upload statistics are kept, and
  **Minimum font size for Shared Files**.

## See also

- [Files](files.md)
- [Downloads and uploads](files-transfers.md)
- [Search](files-search.md)
- [RetroShare links](links.md)
- [Network](network.md)
