# Channels

A channel is a place where its owner publishes posts: a title, a text, a
picture and, most of the time, files to download. Only the owner and the
people who received the publish permissions can post; everybody else reads,
downloads and comments.

Channels spread from friend to friend. You see the channels your friends are
subscribed to, and the channels you subscribe to are offered to your own
friends. Subscribe to a channel to receive its posts.

## Layout

- The title bar shows **Channels** and the **?** help button.
- On the left, the list of channels with its toolbar.
- On the right, the channel that is selected in the list. Several channels can
  be open in tabs.

### The list of channels

The toolbar above the list holds:

- the filter field, which keeps the channels that match the text you type. The
  icon inside the field chooses where to search: **Title** or
  **Description**;
- the **Create Channel** button;
- the **Search entire network...** field: type words and press Enter to look
  for channels beyond the ones your friends know.

The channels are sorted in four categories:

- **My Channels**: the channels you administrate.
- **Subscribed Channels**: the other channels you are subscribed to.
- **Popular Channels**: among the channels you are not subscribed to, the ones
  with the most recent activity.
- **Other Channels**: the remaining channels.

A fifth kind of entry, **Search for "..."**, appears for each network search
and lists its results.

The columns are **Name**, the number of unread posts and the popularity, which
is the number of your friends that are subscribed. Four more columns are
hidden at first: **Friend's Posts** (the number of posts your friends have),
**Last Post**, **Search Score** and **Description**. Hover a channel to see a
summary: popularity, your rights, date of the last post, description and ID.

Click a channel to open it. Click with the middle button to open it in a new
tab.

### An open channel

The toolbar holds the subscribe button, **Add new post**, two view buttons and
a **Search...** field. Below it, three tabs:

- **Channel details**: the logo and name of the channel, its
  **Administrator:**, the date it was **Created:**, the **Last activity:**,
  the **Distribution:** (who can see it), the **Sync period:** and the
  description.
- **Posts**: the posts, as a grid of thumbnails or as a list, and under them
  the selected post in three tabs:
  - **Details**: thumbnail, date, title and text of the post, and a download
    button when the post has files you do not have yet.
  - **Files**: the files of the post with the columns **File**, **Size** and
    **Status**.
  - **Comments**: the comments on the post. The tab title shows their number.
- **All files**: every file of every post of the channel, with one more
  column, **Published**.

Selecting a post marks it as read. In the grid, hold Ctrl and turn the mouse
wheel to change the size of the thumbnails.

## Buttons

### Above the list

- **?** (help button): shows a short introduction to channels.
- **Create Channel**: opens the **Create New Channel** dialog.

### In an open channel

- **Subscribe**: subscribes to the channel. Once subscribed, the button reads
  **Subscribed** followed by the popularity of the channel, and opens a menu
  with **Unsubscribe**. For a channel found by a network search the button
  reads **Channel info missing** until you have requested its data, then
  **Channel info requested...** while it is on its way.
- **Add new post** (tooltip **Post to Channel**): opens the **New Channel
  Post** window. Enabled only when you have the right to publish in the
  channel.
- **Show unread posts only**: when pressed, hides the posts you have read.
  Shown on the **Posts** tab only.
- **Click to switch to list view** / **Click to switch to grid view**:
  switches the posts between thumbnails and a list with one post per line.
  Shown on the **Posts** tab only.
- **Search...**: keeps the posts whose title contains all the words you type,
  and, in **All files**, the files whose name contains them.
- Download button (**Details** tab, tooltip **Download this file:** or
  **Download All these ... files:**): downloads the files of the post that you
  do not have yet. The tooltip lists them with their total size.

### In the Status column of a file

- **Download**: downloads the file to the download directory of the channel.
- While the file downloads, a progress bar is shown with buttons to pause,
  resume and cancel. Cancelling asks for a confirmation and deletes what was
  downloaded.
- **Play** or **Open file**: opens the file once it is complete. The label is
  **Play** for media files.
- **Open folder**: opens the directory that contains the file.

### Comments tab

- **Comment**: opens the **Make Comment** dialog to write a comment on the
  post.
- Sort box: **Hot**, **New** or **Top** orders the comments by activity, date
  or score.
- **Voter ID:** the identity that signs your votes on comments.

## Context menus

### A channel in the list

- **Open in new tab**: opens the channel in a tab of its own.
- **Subscribe** / **Unsubscribe**: starts or stops receiving the posts of the
  channel.
- **Create Channel**: opens the **Create New Channel** dialog.
- **Show Details**: shows the information of the channel in a read-only
  dialog.
- **Edit Details**: opens the **Edit Channel** dialog. Enabled only for the
  channels you administrate.
- **Synchronise posts of last...**: how far back your node asks its friends
  for posts: **5 days**, **2 weeks**, **1 month**, **3 months**, **6 months**,
  **1 year**, **3 years**, **5 years** or **Indefinitly**. The current value
  is greyed out. Enabled only when you are subscribed.
- **Store posts for at most...**: how long posts are kept on your disk, with
  the same choices. Older posts are deleted. Enabled only when you are
  subscribed. The two periods follow each other: RetroShare never
  synchronises further back than it stores.
- **Share publish permissions...**: opens a dialog where you tick the friends
  who will be allowed to post in the channel. Use it also to send the
  permissions to another node of yours. Permissions cannot be taken back.
  Enabled only when you have the publish permissions yourself.
- **Copy RetroShare Link**: copies a link to the channel.
- **Mark all as read**, **Mark all as unread**: changes the state of every
  post of the channel. Enabled only when you are subscribed.
- **Enable Auto-Download** / **Disable Auto-Download**: when enabled, the
  files of new posts are downloaded as they arrive. Only when you are
  subscribed.
- **Set download directory**: where the files of this channel are saved. The
  submenu lists **[Default directory]**, your shared directories and
  **Specify...**, which lets you pick any directory. The current choice is
  marked. Only when you are subscribed.
- **Tree View Options**: **Show Header** shows or hides the header of the
  list, **Show column …** chooses the visible columns, and **Sort by column
  …** is offered when the header is hidden.

### A network search in the list

- **Remove this search**: removes the search and its results.
- **Remove all searches**: removes every search.
- **Request data** (on a result): asks the network for the description of the
  channel. Once it has arrived you can open the channel and subscribe.

### A post

- **Download files**: downloads all the files of the post. Only when the post
  has files.
- **Mark as unread**: only for a post you have read.
- **Copy RetroShare Link**: copies a link to the post.
- **Edit**: opens the post in the **Edit Channel Post** window. Only when you
  have the publish permissions.

### A file

- **Copy Retroshare link**: copies a link to the file. Anybody can paste it in
  RetroShare to download the file.

### A comment

- **Reply to Comment**: opens the **Reply to Comment** dialog.
- **Submit Comment**: writes a new comment on the post.
- **Copy Comment**: copies the text of the comment.
- **Vote Up**, **Vote Down**: gives the comment a positive or negative vote,
  signed by the identity chosen in **Voter ID:**.
- **Show Author**: shows the author of the comment in the [People](people.md)
  page.

### The attachments page of a new post

- **Paste RetroShare Link** / **Paste RetroShare Links**: attaches the files
  whose links are in the clipboard. Greyed out when the clipboard holds no
  file link.

## Dialogs

### Create New Channel

- The picture at the top left is the logo. Click it to choose an image file.
- **Name**: the name of the channel. It is required.
- **Owner:** the identity of yours that is shown as the author of the channel,
  so that readers can contact you.
- **Description**: what the channel is about.
- **Message Distribution**: who can see the channel.
  - **Public**: everybody.
  - **Restricted to Circle**: only the members of the circle you choose in
    the list. See [People](people.md) for circles.
  - **Restricted node group**: only your friend nodes that belong to the
    group you choose.
- **Create**: creates the channel. You are subscribed to it and it appears
  under **My Channels**.

**Edit Details** opens the same dialog, titled **Edit Channel**, with the
button **Update Channel**.

### Channel details

**Show Details** shows, read only: **Name**, **ID**, **Created**, **Posts**,
**Last Post**, **Popularity**, **Author:**, **Distribution:** and the
**Description**.

### Share channel publish permissions

Tick the friends in the **Contacts:** list and confirm. Each of them receives
the right to post in the channel.

### New Channel Post

Two buttons at the top switch between the pages **Channel Post** and
**Attachments**. The second one shows the number of attached files.

On the **Channel Post** page:

- **Channel Post to:** the channel that will receive the post.
- **Add Channel Thumbnail**: chooses the picture of the post. Use the mouse on
  the preview to move and zoom the picture.
- Aspect ratio box: **Auto**, **1:1**, **3:4** or **16:9** sets the shape of
  the thumbnail. **Auto** picks the best fit.
- **Remove Thumbnail** (shown once a picture is set): removes it.
- **Title**: the title of the post. It is required.
- **Hide** / **Show** (arrow button): folds the thumbnail area to leave more
  room for the text.
- The editor below holds the text of the post, with formatting.

On the **Attachments** page:

- **Add File**: attaches files from your disk.
- **Paste retroshare link(s) from clipboard.**: attaches the files whose
  links you copied elsewhere in RetroShare.
- **Remove File**: removes all the attachments.
- **Drag and Drop Files from Search Results**: you can also drop files on the
  window.

**Post** publishes the post. **Cancel** asks whether to discard it.

The files are not stored in the channel: a post only carries their names,
sizes and hashes, and readers download them from the nodes that have them.
Keep the files in one of your shared directories. A file added from outside
your shared directories is shared for 30 days only, and RetroShare warns you
when you post. It also warns when the post refers to files you are not
sharing at all.

**Edit** opens the same window, titled **Edit Channel Post**, with the button
**Update**.

### Make Comment

- **Type in your comment**: the text of the comment. The dialog shows how many
  characters remain.
- **Post as:** the identity that signs the comment. You need an identity to
  comment.
- **Post** publishes the comment.

**Reply to Comment** is the same dialog; it shows **Replying to** and the
author you answer to.

## How to

### Find and follow a channel

1. Look in **Popular Channels** and **Other Channels**, or type words in
   **Search entire network...** and press Enter.
2. For a search result, right click it and choose **Request data**, then wait
   for the channel to arrive.
3. Right click the channel and choose **Subscribe**. The posts arrive from
   your friends over the next minutes.

### Download the files of a post

1. Open the channel and select the post.
2. Click the download button in the **Details** tab, or open the **Files**
   tab and click **Download** next to one file.
3. Follow the transfer in [Downloads and uploads](files-transfers.md).

### Publish a post with a file

1. Put the file in one of your [shared directories](files-shared.md).
2. Open your channel and click **Add new post**.
3. Type a title and a text, and add a thumbnail if you want one.
4. Click **Attachments**, then **Add File**, and choose the file.
5. Click **Post**.

## Related settings

- [Channels](settings-channels.md): tabs, emoticons, maximum size of the
  automatic downloads.
- [Files](settings-files.md): the default download directory.

## See also

- [Forums](forums.md)
- [Boards](boards.md)
- [People](people.md)
- [RetroShare links](links.md)
