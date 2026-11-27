# Downloads and uploads

The **Downloads** tab of the [Files](files.md) page shows the files you are
fetching and, in its lower part, the files that other nodes are fetching from
you. This is where you pause, cancel, reorder and tune your downloads.

## Layout

- At the top, a filter field. Type a text to keep only the downloads that
  match. The magnifying glass of the field chooses what the text is compared
  with: **Name** or **Hash**. Upper and lower case are not distinguished.
- The **Downloads** list, one line per file. The number of downloads is shown
  in the title of the tab.
- Below a separator that you can drag, the **Uploads** list, with the number
  of uploads in its title.

Both lists are updated every second. Click a column header to sort. You can
select several lines with Ctrl or Shift.

### Columns of the Downloads list

- **Name**: the name the file will have once downloaded, with an icon of its
  type.
- **Size**: the size of the complete file.
- **Completed**: the amount already received.
- **Speed**: the current download speed, in KB/s. Empty when nothing is
  coming in.
- **Progress / Availability**: a bar with the percentage received. The bar is
  a map of the file: the parts already received, the parts being downloaded
  and the parts being checked are drawn differently.
- **Sources**: two numbers, such as `2 (5)`: the sources that are sending data
  right now, then in brackets all the sources known for the file.
- **Status**: the state of the download, see below.
- **Speed / Queue position**: for a running download, its priority:
  **Slower**, **Average** or **Faster**. For a queued download, its rank in
  the queue.
- **Remaining**: the amount still to receive.
- **Download time**: the estimated time left at the current speed.
- **Hash**: the hash of the file.
- **Last Time Seen**: the last time data was received for this file.
- **Path**: the directory where the file will be saved. Your incoming
  directory is written `[Download Dir]`.

### Status of a download

- **Waiting**: the download is active but no source is sending data at the
  moment. RetroShare keeps looking for sources.
- **Downloading**: data is being received.
- **Queued**: the download waits for its turn, because the maximum number of
  simultaneous downloads is reached. It starts when a running download ends.
- **Paused**: you paused it. It does not move until you resume it.
- **Checking...**: RetroShare is verifying the data it received against the
  hash of the file. This also happens once at the end of every download,
  before the file is moved to its destination.
- **Complete**: the file is finished and has been moved to its destination
  directory. The line stays in the list until you clear it.
- **Failed**: the download is suspended and cannot go on for the moment.

### Sources

Click the arrow at the left of a download to unfold its sources, one line per
source:

- **Name** shows the name of the friend when the file comes directly from a
  friend. For an anonymous source it shows **Tunnel via** followed by the
  friend through which the tunnel goes. The icon tells friends, anonymous
  tunnels and anonymous end-to-end encrypted tunnels apart.
- **Speed** is the speed of this source alone.
- **Progress / Availability** shows which parts of the file this source has.

The other columns are empty on a source line. An action applied to a source
line applies to its download.

### Columns of the Uploads list

- **Name**: the file that is being sent.
- **Peer**: who receives it: the name of a friend, or **Tunnel via** and the
  name of the friend through which an anonymous tunnel goes. Hover the cell to
  see the tunnel identifier. When several nodes fetch the same file the cell
  shows their number, such as `3 tunnels`, and each of them has its own line
  below the file.
- **Size**: the size of the file.
- **Transferred**: the amount sent to this peer.
- **Speed**: the current upload speed, in KB/s; for a file with several peers,
  the total.
- **Progress**: a bar showing the parts of the file that the peer has got.
- **Hash**: the hash of the file.

## Context menus

### Downloads list

Right click a download, or a selection of downloads. The entries depend on the
status of what is selected.

- **Priority (Speed)...**: shown when at least one selected download is not
  complete.
  - **Priority (Speed)... > Slower**, **Average**, **Faster**: sets how fast
    RetroShare raises the rate it asks from the sources of this file. A
    download on **Faster** reaches a high speed sooner than one on
    **Slower**.
- **Move in Queue...**: shown when at least one selected download is
  **Queued**.
  - **Move in Queue... > Top**: puts the download first in the queue.
  - **Move in Queue... > Up**: moves it one rank up.
  - **Move in Queue... > Down**: moves it one rank down.
  - **Move in Queue... > Bottom**: puts it last.
- **Chunk strategy**: shown when at least one selected download is not
  complete. A file is downloaded in chunks of 1 MB; the strategy is the order
  in which the chunks are asked for.
  - **Chunk strategy > Streaming**: chunks are asked in order, from the
    beginning of the file. Choose it when you want to preview a media file
    while it downloads.
  - **Chunk strategy > Progressive**: the next chunk is taken at random within
    the 50 MB that follow the part already received. A compromise between
    the two others.
  - **Chunk strategy > Random**: chunks are taken at random in the whole file.
    It is the best choice for the other downloaders of the same file, who
    quickly find the parts they miss. On Windows you are warned first, because
    writing in the middle of a large empty file can freeze RetroShare for
    several seconds.
- **Rename file...**: shown for a single download that is not complete. Asks
  for the name the file will have at the end. The characters
  `/ \ | : ? < > *` are not accepted.
- **Set destination directory**: shown when at least one selected download is
  not complete. Chooses where the file will be saved once finished.
  - **Set destination directory > Specify...**: lets you pick any directory.
  - The other entries are your shared directories that exist and are
    writable, under their visible names.
- **Resume**: shown when at least one selected download is **Paused**.
  Restarts it.
- **Pause**: shown when at least one selected download is **Downloading**,
  **Waiting** or **Queued**. Stops it and keeps what was received.
- **Force Check**: shown when no selected download is complete. Marks all the
  chunks already received as to be verified. The chunks that turn out to be
  wrong are downloaded again.
- **Cancel**: shown when no selected download is complete. After a
  confirmation, removes the download and deletes the data already received.
  The Delete key does the same.
- **Play**: shown when a selected file is complete and is a media file. Opens
  it with the default application of your system.
- **Open File**: shown for a single complete download. Opens the file with the
  default application of your system.
- **Preview File**: shown for a single media file that is not complete. Opens
  the part already received in your media player. RetroShare waits until you
  click **OK** in the message **Click OK when program terminates!**, then
  removes the temporary link it created for the player. The preview works
  best with the **Streaming** chunk strategy.
- **Open Folder**: shown for a single download. Opens the directory of the
  file once it is complete, the partials directory before that.
- **Details...**: shown for a single download. Opens the **Details** window
  described below.
- **Clear Completed**: always shown. Removes all the complete downloads from
  the list. The files themselves stay on your disk.
- **Copy RetroShare Link**: shown when at least one download is selected.
  Copies one file link per selected download.
- **Paste RetroShare Link**: shown when the clipboard holds file or directory
  links. Opens the **Download files** window with the files of the links,
  where you choose what to download and where. See [Files](files.md).
- **Expand all**: shown when the list is not empty. Unfolds the sources of
  every download.
- **Collapse all**: shown when the list is not empty. Folds them back.
- **Collection**: always shown. See [Files](files.md) for collections.
  - **Collection > Create Collection...**: creates a collection file from the
    selected downloads.
  - **Collection > Modify Collection...**: enabled when the selection is a
    single complete `.rscollection` file. Opens it in the collection editor.
  - **Collection > View Collection...**: same condition. Opens it in the
    **Download files** window.
  - **Collection > Download from collection file...**: when the selection is a
    single complete `.rscollection` file, opens it in the **Download files**
    window. Otherwise asks for a collection file on your disk.

### Header of the Downloads list

Right click the column headers.

- **Columns**: a submenu with one checkable entry per column, except
  **Name**: **Size**, **Completed**, **Speed**, **Progress / Availability**,
  **Sources**, **Status**, **Speed / Queue position**, **Remaining**,
  **Download time**, **Hash**, **Last Time Seen**, **Path**. Uncheck an entry
  to hide the column. Your choice is remembered.

### Uploads list

- **Open Folder**: shown when a single upload is selected. Opens the directory
  that holds the file.
- **Copy RetroShare Link**: shown when at least one upload is selected. Copies
  one file link per selected upload.
- **Expand all**: shown when the list is not empty. Unfolds the peers of every
  upload.
- **Collapse all**: shown when the list is not empty. Folds them back.

### Header of the Uploads list

- **Columns**: a submenu with the checkable entries **Peer**, **Size**,
  **Transferred**, **Speed**, **Progress** and **Hash**.

## Dialogs

### Details

Opened by **Details...**. The window follows the selection: select another
download and it shows that one. It has two tabs.

**General** shows how the download is going:

- **Chunk map**: one tile per chunk of the file. The legend gives the four
  states: **Done** (received), **Active** (being downloaded),
  **Outstanding** (not received yet) and **Needs checking** (received,
  waiting to be verified).
- **Active chunks**: one bar per chunk being downloaded, with its number and
  how much of it has arrived.
- **Availability map**: one bar per source, showing the chunks that the source
  has, with the number of active sources in the title.
- **File info**: **File name**, **Destination folder**, **File hash**,
  **File size**, **Chunk size**, **Number of chunks**, **Transferred**,
  **Remaining**, **Number of sources**, **Chunk strategy** and
  **Transfer type**. The transfer type is **Anonymous F2F** when RetroShare
  also looks for the file through anonymous tunnels, and
  **Direct friend transfer / Availability assumed** when the sources are
  assumed to have the whole file.

**retroshare link(s)** shows the file link of the download as text.

- **Copy link to clipboard**: copies that text.

## How to

### Watch a video while it downloads

1. Right click the download and choose **Chunk strategy > Streaming**.
2. Wait until the beginning of the file is received.
3. Right click again and choose **Preview File**.
4. When you close the player, click **OK** in the RetroShare message.

### Let one download go before the others

1. If the download is **Queued**, right click it and choose
   **Move in Queue... > Top**.
2. Once it runs, choose **Priority (Speed)... > Faster**.

### Clean the list

1. Right click anywhere in the Downloads list.
2. Choose **Clear Completed**.

## Related settings

- [Files](settings-files.md): **Incoming Directory** and
  **Partials Directory**, **Maximum simultaneous downloads:**,
  **Maximum uploads per friend (0 = no limit)**, **Default chunk strategy:**,
  **Safety disk space limit :**, **Allow direct download:**, and the option
  that downloads the content of a collection file as soon as it arrives.

## See also

- [Files](files.md)
- [Search](files-search.md)
- [Shared files](files-shared.md)
- [RetroShare links](links.md)
- [Statistics](statistics.md)
