# Files

The Files page of the [Preferences](settings.md) sets which directories
RetroShare uses for file sharing, how it scans the directories you share and
how transfers behave.

## Layout

The page has two tabs: **Directories** and **Transfer options**.

## Buttons

### Directories tab, Shared Directories

- **Configure shared directories**: opens the share manager, where you add
  and remove shared directories and choose who can see each of them. See
  [Shared files](files-shared.md).
- **Automatically share incoming directory (Recommended)**: the directory
  where finished downloads are stored is shared as well, so that others can
  get from you the files you have downloaded.
- **Auto-check shared directories every** ... **minute(s)**: RetroShare scans
  the shared directories at this interval to find new and removed files.
  From 1 to 5000 minutes, 10 by default. The scan does not notice a change
  inside an existing file; to force a complete scan, use the check files
  button of the [Shared files](files-shared.md) page. When the box is not
  checked, the directories are no longer scanned automatically and the
  interval is greyed.
- **follow symbolic links**: RetroShare follows the symbolic links it finds
  in the shared directories. Loops and duplicate directories are taken care
  of. When the box is not checked, symbolic links to files and to directories
  are ignored.
- **Ignore duplicate files/directories**: only available when symbolic links
  are followed. A file or directory that is reached a second time through a
  link is silently skipped. This protects the scan against directory loops.
- **Maximum depth (0=unlimited):** the number of directory levels that are
  scanned below each shared directory. Use it to bound the scan when you
  follow symbolic links without ignoring duplicates.
- **Ignore files ending with:** files whose name ends with one of the listed
  suffixes are not shared. Separate the suffixes with `;`, for instance
  `~;.bak;.old`.
- **ignore files starting with:** files whose name starts with one of the
  listed prefixes are not shared. Separate the prefixes with `;`, for
  instance `.;~`.

### Directories tab, Incoming Directory

- The field shows the directory where finished downloads are stored. Click
  the **Browse** button next to it to choose another one. If RetroShare cannot
  use the directory you chose, the previous one is kept and the tooltip of
  the field says **Invalid Input. Have you got the right to write on it?**
- **Automatically donwload RsCollection file content (Not recommended)**:
  when the download of a collection file finishes, RetroShare opens it and
  downloads the files it lists, so you do not review its content first. The
  tooltip warns that a collection may hold a lot of files. Off by default.

### Directories tab, Partials Directory

- The field shows the directory where files are stored while they are being
  downloaded. Click the **Browse** button to choose another one. It cannot be
  a directory that is already shared; in that case the tooltip of the field
  says **Invalid Input. It can't be an already shared directory.**

### Directories tab, Upload Statistics

- **Keep upload statistics for:** in days, from 0 to 3650. Upload statistics
  older than this are deleted. With 0, the default, they are kept forever.
- **Clear upload statistics**: deletes all the upload statistics now.

### Directories tab, bottom

- A text recalls that each shared directory has its own share flags:
  **Browsable by friends** (your friends see the files) and **Anonymously
  shared** (the files can be reached anonymously through distant tunnels).
  You set them in the share manager.
- **Minimum font size for Shared Files**: the size of the text in the
  [Shared files](files-shared.md) page. From 11 to 64; the default is 11 (13
  on macOS).

### Transfer options tab

- **Maximum simultaneous downloads:** the number of downloads that are active
  at the same time. The other downloads wait in the queue. 5 by default.
- **Maximum uploads per friend (0 = no limit)**: the number of files one
  friend can get from you at the same time. 0 by default.
- **Default chunk strategy:** the order in which the parts of a file are
  requested, for new downloads. You can change it for one download in
  [Downloads and uploads](files-transfers.md).
  - **Streaming**: the parts are requested in order from the start of the
    file, which makes it possible to preview the file while it downloads.
  - **Progressive**: the next part is chosen at random within less than
    50 MB after the end of what is already there. A compromise between the
    two others. This is the default.
  - **Random**: the parts are requested in a random order, which is best for
    downloading from several sources. On Windows, RetroShare asks for a
    confirmation first, because writing in the middle of a large empty file
    may hang the program for several seconds.
- **Safety disk space limit :** in MB, from 50 to 1000. When the free space
  of the disk goes below this limit, RetroShare suspends all the transfers
  and stops saving its configuration files, and a window warns you.
- **Max. tunnel req. forwarded per second:** how many tunnel requests of
  other nodes your node forwards per second. Tunnel requests are the small
  messages by which a node looks for an anonymous route to a file. 20 by
  default. With a large bandwidth you may raise it to 30 or 40 so that longer
  tunnels can pass, but these many small packets can slow your own transfers
  down. If you are not sure, keep the default.
- **Allow direct download:** whether a file may be downloaded directly from a
  friend that has it, instead of through anonymous tunnels only.
  - **Yes**: from any friend.
  - **No**: from no friend.
  - **Per user**: only from the friends for whom the option to use them as a
    direct source is checked. See [Friend details](network-details.md).
- **Trust friend nodes with banned files**: RetroShare also takes into
  account the files your friends have banned to decide which files to
  ignore. It is safer to leave it checked.

## How to

### Move your downloads to another disk

1. In the **Directories** tab, click **Browse** in **Incoming Directory**
   and choose the new directory.
2. Do the same in **Partials Directory**, with a directory that is not
   shared.
3. The change applies at once: there is nothing to confirm.

## See also

- [Files](files.md)
- [Shared files](files-shared.md)
- [Downloads and uploads](files-transfers.md)
- [Permissions](settings-permissions.md)
