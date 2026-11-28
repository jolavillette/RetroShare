# Chats

The Chats page lists the chat rooms and shows the room you are talking in. A
chat room is a decentralized group chat: every participant receives every
message, and you take part under one of your identities, not under the name of
your node. You see a room when one of your friend nodes takes part in it, so
you can talk to people who are not your friends.

This page also describes what all chats share: the chat area with its toolbar,
the private chat with a friend node, the distant chat with an identity, and the
history browser. The broadcast chat, which goes to all your connected friends,
is on the [Network](network.md) page.

For chat rooms to work, the clock of your computer must be on time.

## Layout

- On the left, a search field and the **+** button (tooltip **Create new chat
  room**), then the room tree.
- On the right, the room you selected. When you have not entered the selected
  room, or when no room is selected, the right side shows the **Chat Room
  info** pane instead.

### The room tree

The rooms are sorted in four categories:

- **Private Subscribed**: the private rooms you have entered.
- **Public Subscribed**: the public rooms you have entered, with their number.
- **Private**: the private rooms you were invited to and have not entered.
- **Public**: the public rooms your friends take part in and that you have not
  entered, with their number.

The two "Subscribed" categories are hidden when they are empty.

The tree has three columns: **Name**, **Count** (number of participants) and
**Topic**. Only **Name** is shown at first; the two others are switched on in
the context menu. Click a column header to sort.

- A room you have not entered is written in grey.
- A room written in green only accepts identities signed by a profile (no
  anonymous identities).
- The icon of a room changes for a few seconds when someone is typing, enters
  or leaves, and shows a message icon when a new message arrives in a room
  that is not the selected one.
- Hover a room to see its topic, its number of participants, whether auto
  subscribe is enabled, its id and, for a signed room, the line **Security: no
  anonymous IDs**.

Type in the search field (tooltip **Search Chat lobbies**) to keep only the
rooms whose name contains the text.

Double click a room to enter it with your default identity. When you have no
identity yet, the dialog that creates one opens first. When the room requires
signed identities and your default identity is anonymous, the message
**Default identity is anonymous** tells you so and you do not enter: use the
context menu to choose another identity.

### Chat Room info

The pane shows the fields **Chat room Name:**, **Chat room Id:**, **Topic:**,
**Type:** (**Public** or **Private**), **Security:** (**No anonymous IDs** or
**Anonymous IDs accepted**) and **Peers:** (number of participants), followed
by a line that says whether you are subscribed and what a double click does.

### The room view

Once you have entered a room, the right side shows:

- a title bar with your nickname and the room name written as
  `nickname@room`, the search tools, the chat menu button and three buttons
  that belong to the room (see Buttons);
- the conversation;
- the input field and its toolbar;
- on the right, the list of participants with its own search field.

In the list of participants, the bullet in front of a name tells how active the
participant is: green for active during the last 5 minutes, yellow for 5 to 10
minutes, grey for more than 10 minutes, red for a participant you have muted,
blue for yourself. The name of a muted participant is written in red. Hover
the bullet to see for how long the participant has been inactive. Double click
a name to insert `@name` in the input field.

## Buttons

### Chats page

- **?** (help button, top right): opens the manual on this page. The first
  time you open the page, a short introduction to the page is shown in a
  bubble; click the bubble to close it.
- **+** (tooltip **Create new chat room**): opens the **Create Chat Room**
  dialog.

### Title bar of a room

- Dock button (tooltip **Undock to a new window** or **Redock to Main
  window**): moves the room to a window of its own, or brings it back into the
  Chats page.
- Invite button (tooltip **Invite friends to this lobby**): opens the
  **Invite friends** dialog, where you tick the friend nodes that receive an
  invitation to this room.
- Leave button (tooltip **Leave this chat room (Unsubscribe)**): leaves the
  room.

### Title bar of every chat

- Search button (tooltip **Find**, `Ctrl+F`): shows or hides the search field.
  Type at least the number of characters set in the
  [Chats preferences](settings-chat.md), or press `Enter`, to highlight every
  occurrence in the conversation. The tooltip of the search field gives the
  number of items found. Right click the button for the search options
  described in Context menus.
- Arrow buttons (tooltips **Find Previous**, `Ctrl+Shift+G`, and **Find
  Next**, `Ctrl+G`): jump to the previous or next occurrence. They appear
  once a search has found something.
- Mark button (tooltip **Mark this selected text**, `Ctrl+M`): shown while the
  search field is open. Select some text in the conversation and click it to
  remember the place; click it again (tooltip **Return to marked text**) to
  come back to it.
- Notification button, in rooms only: appears when messages of this room are
  waiting in the notification list (for example messages that contain your
  nickname). Its tooltip gives their number. It opens a menu with one entry per
  message, which scrolls to the message, and **Remove All**, which empties the
  list.
- Chat menu button (tooltip **Chat menu**):
  - **Browse Message History**: opens the history browser (see Dialogs).
  - **Save Chat History**: saves the conversation as displayed to a text file.
  - **Clear Chat History**: empties the conversation display. The stored
    history is kept.
  - **Delete Chat History**: asks for confirmation, then empties the display
    and deletes the stored history of this chat.
  - In a room, the last entry is a list of your identities: choose the one
    under which you write in this room. The other participants see a message
    saying that you changed your name. In a room that requires signed
    identities only those are listed.
  - In a private or distant chat, the last entry is **Clear offline
    messages**. It has no effect in this version.

### Information bar

In a private chat, a bar above the conversation tells you when your
correspondent is offline, idle, away or busy. Its **Close** button hides it.
The line under the conversation shows "is typing..." when the other side
writes.

### Toolbar of the input field

- Emoticon button (tooltip **Insert emoticon**): opens the emoticon panel;
  click an emoticon to insert it.
- Sticker button (tooltip **Insert sticker**): opens the sticker panel; click a
  sticker to send it at once as a picture.
- Picture button (tooltip **Attach a Picture**): lets you choose a picture
  file and inserts it in your message. The picture is scaled down so that the
  message stays under the size limit.
- File button (tooltip **Attach a File**): lets you choose one or more files.
  RetroShare hashes them, then inserts a RetroShare link to each file in the
  input field; the files themselves are downloaded by whoever clicks the
  link. A single picture file is inserted as a preview as well. You can also
  drop files on the chat.
- Font button (tooltip **Set font & color**) opens a menu:
  - **Choose font**: sets the font of the selected text, or of all your
    messages in this chat when nothing is selected.
  - **Choose color**: sets the colour of your text in this chat. Font and
    colour are remembered for each chat.
  - **Reset font to default**: goes back to the font set in the
    [Chats preferences](settings-chat.md).
  - **Don't replace tag with Emote Icon.**: when checked, your message is
    sent so that emoticon codes such as `:)` stay as text.
  - **Send as PlainText**: when checked, your message is sent as plain text,
    without font, colour or pictures.
- **Send** button: sends the message. It is disabled when the message is too
  big; its tooltip shows the number of characters left. In a distant chat it
  is disabled while the tunnel is down.

### Keyboard

- `Enter` sends the message; `Ctrl+Enter` or `Shift+Enter` inserts a line
  break. When the option to send with `Ctrl+Enter` is set in the
  [Chats preferences](settings-chat.md), `Ctrl+Enter` sends and `Enter`
  inserts a line break.
- In a room, `Tab` completes the nickname you started to type; press it again
  for the next match, `Shift+Tab` for the previous one. Typing `@`, or
  `Ctrl+Space`, opens the list of participants.
- `Ctrl+F` opens or closes the search, `Ctrl+G` and `Ctrl+Shift+G` go to the
  next and previous occurrence.
- In the conversation, `Ctrl+M` sets or recalls the mark and `Delete` removes
  the selected text from the display.

## Context menus

### Room tree

On a category:

- **Create chat room**: opens the **Create Chat Room** dialog with the
  visibility of the category (public or private) already chosen.

On a room:

- **Leave this room**: only for a room you have entered. Leaves the room and
  switches its auto subscribe off.
- **Enter this chat room**: only for a room you have not entered, when exactly
  one of your identities can be used.
- **Enter this chat room as...**: the same when several identities can be
  used; the submenu lists them with their id.
- **Create an identity and enter this chat room**: shown instead when you have
  no identity. Opens the dialog that creates an identity, then enters the room
  with it.
- **Create a non anonymous identity and enter this room**: shown instead when
  the room requires signed identities and all yours are anonymous.
- **Copy RetroShare Link**: copies a link to the room, which you can paste in
  a chat, a mail or a post.

In a room that requires signed identities, only your signed identities are
counted and listed.

On any item, and on the empty area:

- **Columns > Count**: shows or hides the column with the number of
  participants.
- **Columns > Topic**: shows or hides the topic column.

### List of participants

- **Sort by Activity**: puts the participants who wrote most recently first.
- **Sort by Name**: sorts the participants by name.
- **Start private chat**: opens a distant chat with the participant. Disabled
  on your own identity.
- **Send Message**: opens the mail composer with the selected participants as
  recipients. Disabled on your own identity.
- **Mute participant**: checkable. You no longer see the messages of a muted
  participant in this room. The others still see them. Disabled on your own
  identity.
- **Give positive opinion**, **Give neutral opinion**, **Ban this person (Sets
  negative opinion)**: set your opinion on the identity, which changes its
  reputation (see [People](people.md)). The entry that matches your current
  opinion is disabled.
- **Show author in people tab**: opens the [People](people.md) page on this
  identity. Enabled when exactly one participant is selected.

The entries after the two sort entries are shown only when at least one
participant is selected.

### Conversation

The menu starts with the standard entries to copy text, copy a link location
and select all, then:

- **View Source**: shows the HTML source of the conversation.
- **Save image**, **Copy image**: only on a picture.
- **Clear Chat History**: empties the conversation display.
- **Quote**: only when text is selected. Copies the selection as a quotation
  into the input field.
- **Drop Placemark**: adds a line of dashes at the end of the conversation, to
  mark where you stopped reading.
- **Show Hidden Images**: only on a picture, when pictures are not loaded
  automatically.
- **Import sticker**: only on a picture. Asks for a name and saves the picture
  in your sticker collection.

In a room, a right click on the name of the author of a message adds the
entries of the list of participants, from **Start private chat** to **Show
author in people tab**. Double click the name to insert `@name` in the input
field.

### Input field

After the standard editing entries (undo, redo, cut, copy, paste, delete,
select all):

- **Save image**, **Copy image**: only on a picture you inserted.
- **Paste as plain text**: pastes the clipboard without formatting.
- **Spoiler**: hides the selected text behind a spoiler that the reader opens.
- **Paste RetroShare Link**: pastes the RetroShare links held in the
  clipboard. Disabled when there are none.
- **Paste my certificate link**: inserts a link to your own node certificate.

### Search button

- **Find Case Sensitively**: the search distinguishes upper and lower case.
- **Find Whole Words**: the search only matches whole words.
- **Move To Cursor**: the conversation scrolls to the occurrence found.
- **Don't stop to color after N items found (need more CPU)**: highlights
  every occurrence, however many. N is the limit set in the
  [Chats preferences](settings-chat.md).

### Chat window

Right click an empty part of a private chat window, for example the area
around the avatars:

- **Choose window color...**: lets you choose the background colour of this
  conversation.
- **Give positive opinion**, **Give neutral opinion**, **Ban this person (Sets
  negative opinion)**: in a distant chat only, set your opinion on the identity
  you are talking to.
- **Dock window**: shown twice when chats are grouped in one tabbed window
  (see below). The second entry moves the current tab to a window of its own.
  The first one has no effect in this version.

## Dialogs

### Create Chat Room

- **Room name:** and **Room topic:**: the name and the topic that the other
  users see in their room tree.
- **Identity to use:**: the identity under which you enter the room.
- **Visibility:**: **Public (Visible by friends)** makes the room appear in
  the room tree of your friends, who can enter it and so make it visible to
  their own friends. **Private (Works on invitation only)** shows the room
  only to the friends you invite.
- **Security:** **require PGP-signed identities**: when checked, only
  identities signed by a profile can enter and talk. It prevents anonymous
  spam, since the node of a spammer can then be located by some participants.
  The **Create** button stays disabled while the chosen identity is not a
  signed one.
- **Select the Friends with which you want to group chat.**: tick the friend
  nodes that receive an invitation as soon as the room is created. You can
  leave the list empty and invite later.
- **Create**: creates the room and enters it.

### Invitation to a room

When a friend invites you, the **Join chat room** dialog tells who invites you
and to which room. Choose the identity to use (only signed identities are
offered when the room requires them), then answer **Yes** to enter the room
or **No** to decline. **Yes** is disabled until an identity is chosen.

### Leaving a room

When you close an undocked room window, RetroShare asks **Do you want to
unsubscribe to this chat room?**. Leaving a room removes you from the list of
participants and switches auto subscribe off for this room.

### Message History

Opened by **Browse Message History**. It lists the stored messages of the
chat, with a search field above the list.

- **Copy**: copies the current message to the clipboard as a quotation.
- **Remove**: deletes the selected messages from the stored history.

The context menu of the list has **Mark all** (selects everything), **Copy**,
**Delete** (the selected messages), **Clear history** (the whole stored
history of this chat) and **Quote**, which puts the current message as a
quotation in the input field of the chat and closes the browser. The first
four are disabled when nothing is selected.

Whether a history is kept, and how many messages are shown when a chat opens,
is set for each kind of chat in the [Chats preferences](settings-chat.md).

## Private chat with a friend node

A private chat goes directly to one of your friend nodes, over the encrypted
connection between your two nodes. You start it from the list of friends on
the [Network](network.md) page.

- The chat opens in a window of its own, or as a tab of one shared window
  when **Use a single tabbed window** is checked in the
  [Chats preferences](settings-chat.md).
- The window shows the avatar of your friend and your own. The small button at
  the edge of the avatars (tooltips **Hide Avatar** and **Show Avatar**) folds
  them away. Click your own avatar to change it.
- The title bar shows the name and the status of your friend.
- When your friend is offline, the information bar says so. The messages you
  send are kept and delivered when the friend is online again.

## Distant chat with an identity

A distant chat links one of your identities to someone else's identity through
anonymous tunnels. It is encrypted end to end, and neither side learns the
node of the other. You start it with **Start private chat** in a chat room, or
from the [People](people.md) page.

The window is the one of a private chat, with a coloured light in the title
bar:

- yellow, **Tunnel is being established. Please wait...** or **Tunnel is
  pending**: no tunnel yet. You cannot send.
- green, **End-to-end encrypted conversation established**: you can talk.
- red, **Distant peer has closed the chat**: the other side closed its
  window; or **Connexion refused.**: the other side does not accept distant
  chats.
- grey, **Remote status unknown.**

Closing the window ends the conversation: RetroShare asks **Close
conversation?** and warns that unsent messages are dropped. Who may open a
distant chat with you is set in the [Chats preferences](settings-chat.md).

## How to

### Enter a room and talk

1. Open the **Public** category and double click a room.
2. Type your message in the input field and press `Enter`.

### Create a private room for a few friends

1. Click **+**.
2. Fill in the name and the topic, choose your identity and select **Private
   (Works on invitation only)**.
3. Tick the friends to invite and click **Create**.
4. To invite more friends later, use the invite button in the title bar of the
   room.

### Stop seeing the messages of someone

1. In the list of participants, right click the name.
2. Check **Mute participant**. This applies to this room only.
3. To lower the reputation of the identity everywhere, choose **Ban this
   person (Sets negative opinion)** instead.

### Find a word in a long conversation

1. Press `Ctrl+F` and type the word.
2. Use `Ctrl+G` and `Ctrl+Shift+G` to go through the occurrences.

## Related settings

- [Chats](settings-chat.md): fonts, history, send key, tabbed windows, search,
  distant chat permissions, chat styles.
- [Notify](settings-notify.md): pop-ups and tray notifications for chats and
  rooms.
- [Sound](settings-sound.md): sounds played on new messages.
- [People](settings-people.md): reputation thresholds used for banned
  identities.

## See also

- [Network](network.md): friend list, private chat and broadcast chat.
- [People](people.md): identities, reputation, distant chat.
- [RetroShare links](links.md)
- [Mail](mail.md)
