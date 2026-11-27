# FeedReader

FeedReader is a reader for news feeds. Many web sites (news sites, blogs,
project pages) publish their latest articles as an RSS or Atom feed. You give
FeedReader the address of a feed; it fetches the feed regularly and lists the
new articles, which you read inside RetroShare.

FeedReader can also republish the articles of a feed in a RetroShare forum or
board that you own. Your friends then read the news through RetroShare without
connecting to the web site themselves.

FeedReader is a [plugin](plugins.md): it must be present and enabled in the
[Plugins preferences](settings-plugins.md). It then adds the **FeedReader**
page to the main window. Your node fetches the feeds directly from the web
sites, or through the proxy you set.

## Layout

- On the left, the **Feeds** panel: the tree of your feeds, under the top
  entry **Message Folders**. You can sort feeds into folders. Above the tree
  are two buttons.
- On the right, the messages of the selected feed: a toolbar, the list of
  messages (articles) and, under it, the text of the selected message.
  Several feeds can be open at the same time in tabs.

In the tree, a feed with unread messages is written in bold, followed by the
number of unread messages in brackets. A folder shows the total of its feeds.
While a feed is being updated, its state is shown in brackets after the name:
**Waiting for download**, **Downloading**, **Waiting for process** or
**Processing**. A deactivated feed has a greyed icon, and a feed whose last
update failed has a warning mark on its icon; hover it to read the error. The
button of the FeedReader page in the main window shows the number of new
messages.

Click a feed to show its messages. Click it with the middle button of the
mouse to open it in a new tab. Drag a feed or a folder onto a folder to move
it there.

## Buttons

### Feeds panel

- **Add new feed**: opens the **Create new feed** dialog. The feed is created
  in the selected folder. The button is greyed out when a feed is selected:
  select a folder or **Message Folders** first.
- **Update feed**: fetches the selected feed now. On a folder, it fetches all
  the feeds of the folder.

### Message toolbar

- **Mark messages as read** and **Mark messages as unread**: change the
  state of the selected messages.
- **Mark all messages as read**: marks every message of the feed.
- **Remove messages**: deletes the selected messages. FeedReader remembers
  them so that the next update does not bring them back.
- **Update feed**: fetches this feed now.
- The search field filters the list. Its menu chooses the column searched:
  **Title**, **Date** or **Author**.

### Message list

The columns are **Title**, the read state (a dot), **Date** and **Author**.
Unread messages are in bold.

- Click the dot of a message to switch it between read and unread. The space
  bar does the same for the selected message.
- The Delete key removes the selected messages.

### Message text

- The link button next to the title opens the original article in your web
  browser. Its menu offers **Open link in browser** and **Copy link to
  clipboard**. RetroShare asks for a confirmation before opening a web link,
  see [RetroShare links](links.md).
- The **Hide** / **Expand** button folds and unfolds the text of the message.

A feed that only feeds a forum or a board, without **Local Feed**, keeps no
message in FeedReader. Its list is empty and shows **The messages will be
added to the forum**, **The messages will be added to the board** or **The
messages will be added to the forum and the board**.

## Context menus

### Feeds tree

- **New > Feed**: opens the **Create new feed** dialog, like **Add new
  feed**. Available on a folder and on **Message Folders**.
- **New > Folder**: asks for a name (**Add new folder**, **Please enter a
  name for the folder**) and creates the folder in the selected folder.
- **Open in new tab**: opens the feed in a new tab next to the main one. Not
  shown when the option **Open all feeds in new tab** is on; greyed out on a
  folder and for a feed that is already open in a tab.
- **Edit**: for a feed, opens the **Edit feed** dialog, the same as **Create
  new feed**; for a folder, asks for the new name (**Edit folder**). Not
  available on **Message Folders**.
- **Delete**: removes the feed or the folder after a confirmation (**Do you
  want to remove the feed ...?**). The messages kept for the feed are removed
  with it. Not available on **Message Folders**.
- **Update**: fetches the feed now. Greyed out for a deactivated feed.
- **Deactivate** / **Activate**: stops or resumes the automatic updates of the
  feed. Its messages are kept. Only for feeds, not for folders.

### Message list

- **Mark as read** and **Mark as unread**: for the selected messages.
- **Mark all as read** and **Mark all as unread**: for the whole feed.
- **Copy link**: copies the web addresses of the selected messages to the
  clipboard.
- **Remove**: deletes the selected messages.
- **Add to forum**: with one message selected, lists your forums; choosing
  one opens the forum message editor filled with the title and the text of
  the article, for you to post.
- **Add to board**: the same for your boards, with the board post editor.
- **Retransform**: applies the transformation of the feed again to the
  selected messages. Only for a feed that has a transformation.

## Dialogs

### Create new feed / Edit feed

The **OK** button is available once the dialog holds an address, a name
(unless the name is taken from the feed) and at least one type, with a forum
or a board selected if that type is checked.

Address and name:

- **RSS-Feed-URL:** the web address of the feed.
- **Name:** the name shown in the tree.
- **Description:** a free text.

**Type**: where the articles go. Several boxes can be checked.

- **Local Feed**: the articles are kept in FeedReader, where you read them.
- **Forum**: each article is posted as a message in the forum chosen in the
  list next to the box. The list holds the forums that you administer and
  that have an author identity. The messages are posted under that identity.
- **Update forum information**: keeps the name and description of the forum
  the same as those of the feed.
- **Board**: each article is posted in the board chosen in the list, among the
  boards that you administer.
- **Update board information**: keeps the name and description of the board
  the same as those of the feed.
- **Use first image as board image**: the first picture of the article
  becomes the picture of the board post.
- **Only image**: with the previous option, the board post holds the picture
  alone, without the text of the article, and an article that has no picture
  is not posted.
- **Shrink image**: reduces the picture before posting.

**Misc**:

- **Activated**: the feed is updated automatically. Uncheck it to pause the
  feed.
- **Use name and description from feed**: takes **Name:** and
  **Description:** from what the feed itself announces, and keeps them up to
  date.
- **Embed images**: downloads the pictures of each article and stores them in
  the message. The article can then be read without any connection to the
  web site, which matters for the friends who read it in your forum.
- **Save complete web page (experimental for local feeds)**: instead of the
  summary found in the feed, fetches and keeps the whole web page each article
  links to. Only for a **Local Feed**; while it is checked, **Forum** and
  **Board** cannot be used.

**Authentication (not yet supported)**: the fields **Feed needs
authentication**, **User** and **Password** are greyed out. Feeds that need a
login cannot be read.

**Update interval**:

- **Use standard update interval**: uses the interval of the FeedReader
  preferences.
- **Interval in minutes (0 = manual)**: the interval for this feed. With 0 the
  feed is only fetched when you click **Update**.
- **Last update**: when the feed was last fetched, or **Never**.

**Storage time**:

- **Use standard storage time**: uses the duration of the FeedReader
  preferences.
- **Days (0 = off)**: how long the messages of this feed are kept. Older ones
  are removed. With 0 they are kept for ever.
- **Clear cache**: forgets the messages that you removed from this feed, so
  that an update can bring them back if the web site still lists them.

**Proxy**:

- **Use standard proxy**: uses the proxy of the FeedReader preferences.
- Otherwise, fill in the proxy for this feed: its kind (**HTTP**, **HTTPS**,
  **SOCKS4**, **SOCKS4a**, **SOCKS5** or **SOCKS5h**), its address and its
  port. With **SOCKS4a** and **SOCKS5h** the proxy resolves the name of the
  web site.

**Transformation**:

- **Transformation type** shows the transformation applied to the articles
  of this feed: **No transformation**, **XPath** or **XSLT**.
- **Preview && Transformation**: opens the **Preview** dialog, where the
  transformation is defined. Available once the address is filled in.

### Preview

This dialog fetches the feed with the settings of the dialog you come from and
shows the result, so that you can check a feed before creating it. It is also
where you define a transformation: a rule that keeps or removes parts of each
article, for a feed whose articles come with menus, advertising or other
unwanted blocks.

- **Name:** the name of the feed, followed by the state or the error of the
  download.
- **Previous** and **Next**: go through the messages found; the counter
  between them gives the position.
- **Title:** the title of the displayed message, with its text below.
- **Structure**: shows or hides the tree of the HTML elements of the message,
  which helps you write the expressions.
- **Transformation type**:
    - **No transformation**: the articles are kept as they are;
    - **XPath**: two lists of XPath expressions. **XPath use** keeps only the
      parts of the article that match one of its expressions. **XPath
      remove** deletes the parts that match. Right click a list for **Add**,
      **Edit** and **Delete**;
    - **XSLT**: an XSLT style sheet, typed in the text field, rewrites the
      article. It is applied when you leave the field or press Ctrl+Enter.
- The result of the transformation is shown at once. An error such as **Wrong
  XPath expression**, **Empty XPath result** or **XSLT format error** is
  displayed above the message.

Confirm the dialog to keep the transformation in the feed dialog.

## The FeedReader preferences page

The plugin adds a **FeedReader** page to the preferences, with the values
used by every feed that does not have its own.

- **Update** > **Interval in minutes (0 = manual)**: how often the feeds are
  fetched.
- **Storage time** > **Days (0 = off)**: how long the messages are kept.
- **Proxy** > **Use proxy**: fetches the feeds through the proxy filled in
  next to it (kind, address, port). On a hidden node, use the proxy of Tor or
  I2P with the **SOCKS5h** kind, where the proxy resolves the names of the web
  sites: with plain **SOCKS5** your computer looks the names up itself and
  reveals which sites you read.
- **Open all feeds in new tab**: a click on a feed opens it in its own tab
  instead of replacing the content of the main tab.
- **Set message to read on activate**: a message is marked as read as soon as
  you select it.
- **Save configuration in background (for slow systems, more memory
  needed)**: writes the feeds and their messages to disk without blocking
  RetroShare while it does so.

## Notifications

- The button of the **FeedReader** page shows the number of new messages, and
  the same count can be shown in the system tray like for the other pages:
  see the entry **FeedReader Message** in the
  [Notify preferences](settings-notify.md), tab **Systray**.
- **Feed Reader**, in the list of the **Activity** tab of the
  [Notify preferences](settings-notify.md): when checked, each new article is
  added to the [Activity](activity.md) page. The entry shows the feed name,
  the date, the title and the text, with three buttons: **Expand** shows the
  text, **Set as read and remove item** marks the article as read and takes
  it off the page, **Remove Item** only takes it off the page. The link button
  opens the article in the browser or copies its address.

## How to

### Follow a web site

1. Find the address of the feed on the web site (often behind an "RSS" or
   "Atom" link).
2. In the FeedReader page, select **Message Folders** and click **Add new
   feed**.
3. Paste the address in **RSS-Feed-URL:**.
4. Keep **Use name and description from feed** and **Local Feed** checked and
   click **OK**.
5. Select the feed and click **Update feed** to fetch it at once.

### Publish a feed in one of your forums

1. Create the forum in the [Forums](forums.md) page, with one of your
   identities as author.
2. Create or edit the feed and check **Forum**, then choose the forum in the
   list.
3. Check **Embed images** if the articles should keep their pictures.
4. Uncheck **Local Feed** if you do not want to keep a copy in FeedReader.
5. Click **OK**. The new articles are posted in the forum at each update.

## Related settings

- [Plugins](settings-plugins.md): enabling the FeedReader plugin.
- [Notify](settings-notify.md): Activity entries and tray notifications.

## See also

- [Plugins](plugins.md)
- [Forums](forums.md)
- [Boards](boards.md)
- [Activity](activity.md)
