# Forums

A forum is a place for discussions. Anybody who is subscribed can start a
thread and reply to a message, signing with one of their identities. Forums
spread from friend to friend, like [channels](channels.md): you see the forums
your friends are subscribed to, and your own subscriptions are offered to your
friends.

## Layout

- The title bar shows **Forums** and the **?** help button.
- On the left, the list of forums with its toolbar.
- On the right, the forum that is selected in the list. Several forums can be
  open in tabs.

### The list of forums

The list works like the list of channels; see
[Channels](channels.md) for the columns, the filter field with its **Title**
and **Description** options, and the tooltips. The differences are:

- the categories are **My Forums**, **Subscribed Forums**, **Popular Forums**
  and **Other Forums**;
- the button of the toolbar is **Create Forum**;
- there is no field to search the whole network;
- the number of unread messages counts the replies as well as the threads.

### An open forum

From top to bottom:

- A toolbar with the **Subscribe** button (only while you are not
  subscribed), **New thread**, the name of the forum, three view buttons and
  the **Search forums** field.
- The list of messages, with the columns **Title**, **UnRead**, **Date**,
  the distribution icon and **Author**.
- A second toolbar with **Reply**, the navigation buttons, the date of the
  selected message, the list of its versions when it has been edited, and
  **By** followed by its author.
- The text of the selected message.

When no message is selected, the lower part shows the description of the
forum: name, description, number of subscribers, number of posts at your
friends, last post, synchronisation and storage periods, distribution, owner,
anti-spam rule and moderators. Click the name of the forum in the toolbar to
come back to this description.

In the list of messages:

- Unread messages are shown in bold, and so is a message that has unread
  replies.
- Pinned threads stay at the top of the list and are shown in a different
  colour.
- The distribution icon tells whether the message is forwarded to your
  friends. Hover it for the reason: the author is banned, has no positive
  reputation in a forum that requires one, or the information about the
  identity is missing.
- The title of a message written by an identity you banned is replaced by
  **[ ... Redacted message ... ]**.
- Click a column header to sort.

## Buttons

- **?** (help button): shows a short introduction to forums.
- **Create Forum** (above the list): opens the **Create New Forum** dialog.
- **Subscribe**: subscribes to the forum. Your node then collects the
  messages from your friends and offers the forum to them. To unsubscribe, use
  the context menu of the forum in the list.
- **New thread** (tooltip **Start new Thread for Selected Forum**): opens the
  message window to start a new discussion. Enabled only when you are
  subscribed.
- **Threaded**: shows the replies under the message they answer, as a tree.
- **Flat**: shows all the messages at the same level.
- **Latest post in thread**: shows the tree and sorts the threads by the date
  of their most recent reply instead of the date of their first message.
- **Search forums**: keeps the messages that match the text. The icon inside
  the field chooses the column that is searched: **Title**, **Date** or
  **Author**. The tooltip of the field gives the number of results.
- **Reply** (tooltip **Reply Message**): opens the message window to answer
  the selected message. Enabled when you are subscribed and a message is
  selected.
- **Previous Thread**, **Next Thread**: selects the previous or next message
  at the same level.
- **Next unread message**: jumps to the next unread message, opening the
  threads on the way.
- **Download all files**: starts the download of every file link found in the
  text of the message. Enabled only when the message contains links.
- Version list (shown only for a message that has been edited): **(Latest)**
  is the current text; choose an **(Old)** entry to read a previous version.

Click in the **UnRead** column of a message to switch it between read and
unread.

## Context menus

### A forum in the list

The entries are the same as for a channel and are described in
[Channels](channels.md):

- **Open in new tab**
- **Subscribe** / **Unsubscribe**
- **Create Forum**
- **Show Details**
- **Edit Details**: enabled only for the forums you administrate.
- **Synchronise posts of last...** and **Store posts for at most...**: enabled
  only when you are subscribed.
- **Copy RetroShare Link**
- **Mark all as read**, **Mark all as unread**: enabled only when you are
  subscribed.
- **Tree View Options**: **Show Header**, **Show column …** and, when the
  header is hidden, **Sort by column …**.

A forum has no **Share publish permissions...** entry and no download
entries: everybody who is subscribed can post, and a forum carries no files
of its own.

### A message

- **Edit**: opens the message in the message window. Shown for your own
  messages, and for any message when one of your identities is a moderator of
  the forum. Not shown for a pinned message. The edited message replaces the
  original; the previous versions stay readable from the version list.
- **Pin this post up** / **Un-pin this post**: keeps the thread at the top of
  the list for every reader, or releases it. Only for the administrator of the
  forum, on the first message of a thread.
- **Reply**: answers the message. Greyed out when you are not subscribed.
- **Start New Thread**: starts a new discussion. Greyed out when you are not
  subscribed.
- **Copy RetroShare Link**: copies a link to the message.
- **Mark as read** / **Mark as unread**: only the one that applies is shown.
- **Mark as read (with children)**, **Mark as unread (with children)**: do the
  same for the message and all its replies. Shown for a message that has
  replies.
- **Expand all**, **Collapse all**: opens or closes every thread.
- **Expand subtree**: opens all the replies under the message. Shown for a
  message that has replies.
- **Author's reputation > Give positive opinion**, **Give neutral opinion**,
  **Give negative opinion**: sets your opinion of the author. A negative
  opinion bans the author: the messages are hidden and no longer forwarded to
  your friends. The opinion you already have is not listed. Not shown for
  your own identities. See [People](people.md).
- **Show author in people tab**: shows the author in the
  [People](people.md) page.
- **Reply to author with private message**: opens a new
  [mail](mail.md) addressed to the author, quoting the message. Greyed out
  when you are not subscribed.

### The header of the message list

- **Title**, **Read**, **Date**, **Distribution**: show or hide the column.
- **Show text from banned persons**: when ticked, the text of messages
  written by identities you banned is displayed anyway.

## Dialogs

### Create New Forum

- **Name**: the name of the forum. It is required.
- **Owner:** the identity of yours that is published as the contact of the
  forum and has the rights of a moderator. You can leave it empty to
  administrate the forum anonymously.
- **Add moderators**: shows the **Moderators:** list, where you tick the
  identities that can edit and pin the messages of other people. The box
  above the list restricts it to **All People** or **My Contacts**.
- **Description**: what the forum is about.
- **Message Distribution**: who can see the forum.
  - **Public**: everybody.
  - **Restricted to Circle**: only the members of the circle you choose. See
    [People](people.md) for circles.
  - **Restricted node group**: only your friend nodes that belong to the
    group you choose.
- **Spam-protection**, **Posts permissions:** decides how messages from
  unknown authors travel.
  - **All allowed**: every message is forwarded.
  - **Defavor unsigned IDs**: a message signed by a pseudonym is forwarded
    only when its author has a positive reputation.
  - **Defavor unsigned IDs and IDs from unknown nodes**: the same rule also
    applies to identities signed by a profile your node does not know.
- **Create**: creates the forum. It appears under **My Forums**.

**Edit Details** opens the same dialog, titled **Edit Forum**, with the button
**Update Forum**.

### Forum details

**Show Details** shows, read only: **Name**, **ID**, **Created**, **Posts**,
**Last Post**, **Popularity**, **Author:**, **Distribution:**,
**Moderators:**, **Anti Spam:** and the **Description**.

### Start New Thread, Post Forum Message, Edit Message

The message window has the same fields for a new thread, a reply and an
edit.

- **Forum**: the forum that receives the message. For a reply, the title of
  the message you answer is added.
- **Post as**: the identity that signs the message. Only the identities
  allowed to post are offered; RetroShare tells you when none of yours
  qualifies, because of the circle of the forum or of its anti-spam rule.
- **Title**: the title of the message. A reply starts with the title of the
  message it answers.
- **Text**: the message. A counter warns you when the message is too big.
- **Attach File**: adds files from your disk. They are hashed and a link to
  each of them is inserted in the text; readers download them from your
  node, so keep them shared.
- The emoticon button inserts an emoticon.
- **Attach a Picture**: inserts a picture in the text.
- **Post** publishes the message; the button reads **Update** for an edit.
  **Cancel** asks for a confirmation when the message is not empty.

You can also drop files on the window to attach them.

## How to

### Take part in a discussion

1. Select the forum in the list and click **Subscribe**.
2. Wait for the messages to arrive from your friends.
3. Select a message and click **Reply**, or click **New thread**.
4. Choose the identity in **Post as**, write the message and click **Post**.

### Keep spam out of a forum you read

1. Right click a message of the spammer.
2. Choose **Author's reputation > Give negative opinion**.

The author's messages disappear from every forum, channel and board, and your
node stops forwarding them.

### Create a private forum

1. Create a circle with the people you want; see [People](people.md).
2. Click **Create Forum**, choose **Restricted to Circle** and select the
   circle.
3. Click **Create**. Only the members of the circle receive the forum.

## Related settings

- [Forums](settings-forums.md): tabs, emoticons, embedded images, marking
  messages as read, expanding new messages.
- [People](settings-people.md): reputation rules.

## See also

- [Channels](channels.md)
- [Boards](boards.md)
- [People](people.md)
- [RetroShare links](links.md)
