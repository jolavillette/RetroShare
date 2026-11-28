# RetroShare links

A RetroShare link is a piece of text that starts with `retroshare://` and
points to something inside the RetroShare network: a file, a forum, a channel
post, a person, a chat room. You copy a link in one place and paste it in a
chat, a mail or a forum post. Whoever reads the message clicks the link and
RetroShare takes them to the subject.

A link looks like this:

`retroshare://file?name=Holiday.jpg&size=482113&hash=91a6...`

In a message it is displayed as a clickable name, like a link on a web page.
A link contains no address of your computer. It only works inside RetroShare:
a web browser cannot open it, unless RetroShare is registered as the program
for these links (see below).

## Kinds of links

### Files

- **File**: the name, the size and the hash (fingerprint) of a file. Clicking
  it opens the **Download files** dialog, where you confirm what to download
  and where to. The file is then looked for among your friends and through
  anonymous tunnels. If you already have the file and you clicked a single
  link, RetroShare offers to open it instead: it warns you first (**Warning:
  Retroshare is about to ask your system to open this file.**) and opens the
  file only if you answer **Yes**. A collection file (`.rscollection`) that
  you already have is opened in the collection dialog.
- **Folder** (`retroshare://collection?...`): a whole directory with its
  sub-directories, with the number of files and the total size. Clicking it
  opens the **Download files** dialog with the tree of the directory, where
  you choose the files you want.
- **File from one friend** (`retroshare://extra?...`): a file link that also
  names the node that offers the file. It starts the download from that node.

Characters that are not allowed in a file name (`/ \ " * : ? < > |`) are
replaced by `_` and a warning tells you so.

### People

- **Certificate** (`retroshare://certificate?...`): the RetroShare ID of a
  node. Clicking it opens the friend wizard with this ID, so that you can add
  the node as a friend. See [Adding friends](start-friends.md).
- **Person** (`retroshare://person?...`): a profile, by its PGP key
  identifier. Clicking it opens the details of the profile if your node
  already knows it; otherwise the result window says **Friend not found**.
- **Identity** (`retroshare://identity?...`): an identity (pseudonym) with its
  public data. Clicking it adds the identity to your [People](people.md) page
  and confirms with **Identity added to People**. You can then chat with it
  or write to it.
- **Circle** (`retroshare://circle?...`): opens the [People](people.md) page
  on this circle.
- **Message** (`retroshare://message?...`): opens the mail composer with the
  person as recipient, and a subject if the link carries one. The recipient is
  a friend node or an identity. It does nothing but warn you when the identity
  is one of yours.

### Content

- **Forum** (`retroshare://forum?...`): opens the [Forums](forums.md) page on
  the forum, or on one message of the forum when the link names a message.
- **Channel** (`retroshare://channel?...`): opens the
  [Channels](channels.md) page on the channel, or on one post.
- **Board** (`retroshare://posted?...`): opens the [Boards](boards.md) page on
  the board, or on one post.
- **Chat room** (`retroshare://chat_room?...`): opens the [Chats](chats.md)
  page on the room. A link to the group chat opens the
  [Network](network.md) page.
- **Search** (`retroshare://search?keywords=...`): opens the
  [Search](files-search.md) tab and starts a search for the keywords.

A forum, channel, board or chat room link only works if your node already
knows the forum, channel, board or room. Otherwise the result window says, for
example, **Forum not found** or **Channel message not found**: the content has
not reached you yet through your friends.

### After a click

When something could not be done, a window titled **Result** sums up what
happened to the links: how many were processed, and for the others the reason
(**File exist**, **Friend not found**, **Forum not found**, **Malformed
links**, **Invalid links**...).

## Creating a link

Almost every list of RetroShare has an entry to copy the link of the selected
item to the clipboard. The most common ones:

- **Copy RetroShare Link**: in the context menu of a download, of a search
  result, of your shared files and directories, of a forum, a forum message, a
  channel, a channel post, a board, a board post, a chat room, and of a
  search in the list of your searches.
- **Send RetroShare Link**: in search results and shared files, opens a new
  mail that already contains the link.
- **Copy certificate link**: in the friend list of the [Network](network.md)
  page, copies the link to the ID of a friend node.
- **Copy identity to clipboard** and **Copy retroshare link**: in the
  [People](people.md) page, for an identity and for a circle.

Several items selected together give several links, one per line.

When you attach a file in a chat window or in a forum post, RetroShare
computes its hash and inserts the file link for you.

## Pasting a link

In every field where you write a message (chat, mail, forum post, channel
post, comment), right click to get:

- **Paste RetroShare Link**: inserts the links that are in the clipboard as
  clickable links. The entry is greyed out when the clipboard holds no link.
- **Paste my certificate link**: inserts the link to your own RetroShare ID,
  a quick way to give it to someone.

Other places accept links directly:

- **Paste RetroShare Link** in the context menu of the downloads list (shown
  when the clipboard holds a file or folder link): starts the download of the
  files, as if you had clicked the links. See
  [Downloads and uploads](files-transfers.md).
- **Paste certificate link** in the friend list of the
  [Network](network.md) page: opens the friend wizard with the ID in the
  clipboard.
- **Paste RetroShare Link** or **Paste RetroShare Links** in the list of
  attached files when you write a channel post: attaches the files of the
  links to the post.
- **Paste RetroShare Link** in the list of recommended files of the mail
  composer.

## Web links

RetroShare also handles the ordinary `http://` and `https://` links found in
messages. Before handing one to your web browser, it asks **Do you want this
link to be handled by your system?** and shows the full address: check that
the link has not been forged to lead you to a malicious site. The dialog has
two check boxes:

- **Don't ask me again**: no more question until RetroShare is closed.
- **Never ask me again**: no more question at all.

## Opening links from outside RetroShare

### Registering the retroshare:// protocol

When the protocol is registered, a click on a `retroshare://` link in another
program (a web browser, a mail program) is sent to RetroShare.

- On Windows, check **Register retroshare:// as URL protocol** in the
  [General preferences](settings-general.md). The option is available when
  you have the rights to change the registry; the note next to it tells you.
- On Linux, the protocol is declared by the desktop file installed with the
  RetroShare package. Nothing is to be done in RetroShare. If the file is
  missing or does not point to the program you are running, the General
  preferences show **!!!The RetroShare's desktop file is missing or wrong!!!**
  and you should install RetroShare with a package.

### From the command line

A link or a collection file given to the `retroshare-gui` command is passed to
the RetroShare that is already running, which opens it as if you had clicked
it:

`retroshare-gui "retroshare://file?name=...&size=...&hash=..."`

This needs **Accept operating systems calls** to be checked in the
[General preferences](settings-general.md) of the running RetroShare (it is
checked by default). See [Command line options](commandline.md).

## Related settings

- [General](settings-general.md): **Register retroshare:// as URL protocol**,
  **Accept operating systems calls**.

## See also

- [Downloads and uploads](files-transfers.md)
- [Shared files](files-shared.md)
- [Search](files-search.md)
- [Command line options](commandline.md)
