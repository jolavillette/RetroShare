# Search

The **Search** tab of the [Files](files.md) page finds files by name: in your
own shared files, in the directories your friends let you browse, and in the
rest of the network through anonymous searches. From the results you start a
download, copy a link or mark a file as bad.

## Layout

- At the top, the search field, the **Search** and **Advanced** buttons, the
  three check boxes that set where to look, the **Clear** button and the
  banned files button.
- On the left, the list of the searches you have made, with the columns
  **KeyWords** and **Results** (the number of files found so far). Click a
  search to see its results. A search that was started with a file type shows
  the type in brackets after the keywords.
- On the right, the results of the selected search.
- Under the results, the **Filter:** field, **Max results:**, the file type
  list and the **Download selected** button.

Results keep arriving for a while after you start a distant search: the
answers come back from nodes that are several hops away.

### Columns of the results

- **Filename**: the name of the file, with an icon of its type. A directory
  found in the files of a friend can be unfolded to see what it contains.
  When several sources give different names to the same file, the other
  names are listed below it.
- **Size**: the size of the file.
- **Sources**: two numbers separated by a slash, such as `1/4`: the number of
  friends that have the file, then the number of answers received from
  anonymous searches. Hover the cell to see through which friend the first
  answer came.
- **Type**: the kind of file, from its extension. **Folder** for a directory.
- **Age**: the time since the file was last modified, as known by its source.
- **Hash**: the hash of the file. Empty for a directory.

The colour of a line tells you about the file: files you already have, files
you are downloading, and files with few or many sources each have their own
colour, which depends on the style of the interface.

## Buttons

- Search field (tooltip **Enter a keyword here (at least 3 char long)**):
  type one or several words. With several words, a file matches when its name
  contains all of them. Upper and lower case are not distinguished. Press
  Enter to start the search.
- **Search** (tooltip **Start Search**): starts the search. Disabled until the
  text is at least 3 characters long.
- **Advanced** (tooltip **Advanced Search**): opens the advanced search
  window described below.
- **Browsable files** (tooltip **Search inside "browsable" files of your
  friends**): looks in the directories that your friends let you see in the
  **Friends files** tab. The answer is immediate, because the lists of these
  files are already on your node.
- **Distant** (tooltip **Multi-hop search at distance 6 in the network (always
  reports available files)**): sends an anonymous search request through your
  friends, to their friends and so on, up to 6 hops away. Only files that
  their owners allow to be searched anonymously are found.
- **Own files** (tooltip **Include files from your own file list in the search
  result**): also looks in the files you share.
- **Clear** (tooltip **Close all Search Results**): removes all the searches
  and their results.
- Banned files button (red biohazard icon, no label): opens the list of banned
  files described below.
- **Filter:** (tooltip **Filter Search Result**): type a text to keep only the
  results whose file name contains it.
- **Max results:**: the largest number of different files kept for one
  search, from 100 to 5000. Answers beyond that number are ignored. The value
  is remembered.
- File type list: shows only the results of one kind. The choices are **Any**,
  **Archive**, **Audio**, **CD-Image**, **Document**, **Picture**,
  **Program**, **Video** and **Directory**. The kind is deduced from the
  extension of the file name. The type selected when you start a search is
  remembered with that search.
- **Download selected** (tooltip **Download Selected**): downloads the
  selected results. A double click on a result does the same.

## Context menus

### Result

Right click one or several selected results.

- **Download**: starts the download of the selected files. For a directory,
  all the files it contains are downloaded into a directory of the same name
  inside your incoming directory. Files that you already have are skipped, and
  the message **Skipping Local Files** tells you so. The downloads appear in
  the [Downloads](files-transfers.md) tab.
- **Copy RetroShare Link**: copies one file link per selected file.
  Directories are left out.
- **Send RetroShare Link**: opens a new mail that contains the links of the
  selected files. See [Mail](mail.md).
- **Mark as bad**: adds the selected files to your list of banned files and
  removes them from the results. See the banned files list below.
- **Open Folder**: shown when a single file that you already have is
  selected. Opens the directory that holds it.
- **Collection**: see [Files](files.md) for collections.
  - **Collection > Create Collection...**: creates a collection file from the
    selected files.
  - **Collection > Modify Collection...**: enabled when the selection is a
    single `.rscollection` file that you already have.
  - **Collection > View Collection...**: same condition.
  - **Collection > Download from collection file...**: opens a collection
    file to download its content.

### Search (left list)

Right click a search.

- **Search again**: removes the search and starts it again with the same
  keywords and the same file type. Disabled for an advanced search.
- **Remove**: removes the search and its results. The Delete key does the
  same.
- **Remove All**: removes all the searches, like the **Clear** button.
- **Copy RetroShare Link**: copies a search link that holds the keywords.
  Whoever clicks it starts the same search on their own node. Disabled for an
  advanced search. See [RetroShare links](links.md).

## Dialogs

### RetroShare: Advanced Search

The advanced search combines several conditions on the name, the location,
the date or the size of the files. It looks in the places checked in the
Search tab: **Browsable files**, **Distant**, **Own files**. The window closes
when you leave the Search tab.

The **Search Criteria** box holds one line per condition. A line is made of:

- on every line but the first, how the line is combined with what is above:
  - **and**: both must be true;
  - **and / or**: at least one must be true;
  - **or**: one must be true, but not both;
- what is tested: **Name** (the file name), **Path** (the directories that
  lead to the file), **Extension**, **Hash**, **Date** (the date of last
  modification) or **Size**;
- the comparison:
  - for **Name**, **Path**, **Extension** and **Hash**: **contains** (at least
    one of the words you type), **contains all** (every word) or **is** (the
    exact text);
  - for **Date** and **Size**: **less than**, **less than or equal**,
    **equals**, **greater than or equal**, **greater than** or
    **is in range**;
- the value:
  - a text, with the **ignore case** check box, checked by default. For
    **Hash** only the characters of a hash (0-9, a-f) are accepted and there
    is no check box;
  - a date, written year-month-day;
  - a number and its unit, **KB**, **MB** or **GB**;
  - with **is in range**, two dates or two sizes separated by **to**;
- on every line but the first, a button that removes the line (tooltip
  **Delete this expression**).

The buttons of the window:

- Add button (tooltip **Add a further search criterion.**): adds a line.
- Reset button (tooltip **Reset the search criteria.**): removes all the lines
  and starts again with a single empty one.
- **Cancel** (tooltip **Cancels the search.**): closes the window.
- **Search** (tooltip **Perform the advanced search.**): closes the window and
  starts the search. A new line appears in the list of searches with a summary
  of the conditions.

### Banned files

Opened by the banned files button. A banned file is a file that you decided
to keep out of your part of the network: your node no longer forwards search
results nor data for it to your friends. The list of your banned files is
shared with your friends, and their nodes apply it too unless they have
turned the option off. A file disappears from the whole network only if a
large part of the users ban it.

The list has the columns **Filename**, **Hash**, **Size** and
**Banned since...**.

- **Remove** (right click a file): takes the file out of the list.

## How to

### Find and download a file

1. Check **Distant**, and **Browsable files** if you also want the files of
   your friends.
2. Type at least 3 characters and press Enter.
3. Wait for the results, then sort them by **Sources**: a file with more
   sources downloads better.
4. Double click the file.

### Find large recent videos

1. Click **Advanced**.
2. On the first line choose **Extension**, **contains**, and type the
   extensions, for instance `mkv mp4`.
3. Click the add button. On the new line choose **and**, **Size**,
   **greater than**, type `500` and choose **MB**.
4. Click **Search**.

### Ban a file

1. Right click the file in the results.
2. Choose **Mark as bad**.

## Related settings

- [Files](settings-files.md): **Trust friend nodes with banned files**
  decides whether your node also applies the lists of banned files of your
  friends.

## See also

- [Files](files.md)
- [Downloads and uploads](files-transfers.md)
- [Shared files](files-shared.md)
- [RetroShare links](links.md)
