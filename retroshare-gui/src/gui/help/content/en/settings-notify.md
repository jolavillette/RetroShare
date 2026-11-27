# Notify

RetroShare can tell you about what happens on your node in four ways: an item
in the [Activity](activity.md) page, a toaster (a small window that pops up in
a corner of the screen for a few seconds), an icon in the system tray, and a
message balloon of the system tray. This page of the preferences chooses, event
by event, which of these are used.

Every change is applied at once.

## Layout

The page has three tabs: **Activity**, **Toasters** and **Systray**.

## Buttons

### Activity tab, group News Feed

Each checkbox enables one family of items in the [Activity](activity.md) page.
An unchecked family creates no item; the items already in the page stay.

- **Friend Connected**: an item when a friend node connects. The same box also
  enables the item that reports a friend whose clock differs from yours, and
  the item shown when a friend sees your node under another IP address than
  the one you announce. Off by default.
- **Circles**: membership requests, invitations, members who join or leave, and
  revoked memberships in the circles that concern you. On by default.
- **Channels**: new or updated channels, new or updated posts, new comments,
  and the reception of the publish rights of a channel. On by default.
- **Forums**: new or updated forums, new or updated messages, and changes in
  the list of moderators of a forum. On by default.
- **Boards**: new or updated boards, new posts and new comments. On by default.
- **Wire**: activity of the experimental service The Wire. The box has an
  effect only in a build that includes this service.
- **Mails**: an item for each new mail. On by default.
- **Chats**: off by default. No item of the current Activity page depends on
  this box.
- **Friend Requests**: an item when a connection attempt fails the security
  checks: a node that is not your friend tries to connect (a friend request),
  a certificate is missing, wrong or badly signed, or the address is
  blacklisted. On by default.
- **Ip security**: an item titled **Connection refused by peer** when a friend
  node refuses your connection. On by default.
- **Non-Retroshare connections**: one item that collects the addresses of the
  connections received by your port that do not identify themselves as a
  RetroShare node. On by default.

The last three boxes work together: security events are looked at only when
**Friend Requests** or **Non-Retroshare connections** is checked. To silence
them all, uncheck both.

A plugin that can feed the Activity page adds its own checkbox under these.

### Activity tab, group Message

- **Friend Requests**: when checked, a failed connection attempt that creates
  a security item, such as a friend request from a node that is not your
  friend, also puts a message titled "... wants to be friends with you on
  RetroShare" in the inbox of your [Mail](mail.md) page. A request that is
  already in the inbox is not added twice. On by default.

### Toasters tab

- **Disable All Toasters** (tooltip **Disable All Toaster temporarily**): a
  button that stays pressed. While it is pressed no toaster is shown, whatever
  the boxes below say. The choice is not saved: toasters are enabled again the
  next time RetroShare starts. The toaster button of the status bar does the
  same thing.

Group **Toasters**: each box enables the toaster of one event.

- **Friend Connect**: a friend node connects. On by default.
- **New Message**: a new mail arrives. On by default.
- **Download completed**: a download finishes.
- **Private Chat**: a message arrives in a private chat. No toaster is shown
  when the chat window of that conversation is the active one.
- **Group Chat**: a message arrives in the broadcast chat with all your
  friends. No toaster is shown when you are looking at that chat.
- **Chat Room**: a message arrives in a chat room. No toaster is shown when the
  room is the active chat, nor for a participant you have muted.
- **Friend Requests**: a node that is not your friend tries to connect.

A plugin that has toasters adds its own boxes under these, alone or as a
titled group of boxes.

Group **Toaster position**:

- **Position**: the corner of the screen where toasters appear: **Top Left**,
  **Top Right**, **Bottom Left** or **Bottom Right** (the default).
- **X Margin** and **Y Margin**: the horizontal and vertical distance, in
  pixels, between the toaster and the edge of the screen. From 0 to 100, 10 by
  default.
- **Test**: shows a sample toaster for every box that is checked, at the
  position and with the margins currently set. Use it to check the placement
  on your screen.

### Systray tab, group Systray Icon

This group has one line per kind of unread element: **Message** (mail),
**Download completed**, **Private Chat**, **Chats** (chat rooms),
**Channel Post**, **Forum Post**, **Board Post**, and the lines added by
plugins or by optional services. Each line has three boxes:

- the box with the name: when checked, an icon appears in the system tray as
  long as there are unread elements of this kind, with their count in its
  tooltip. Clicking the icon opens the matching page. On by default. When it is
  unchecked the two other boxes are cleared and greyed.
- **Combined**: instead of its own tray icon, this kind of element is reported
  by the main RetroShare tray icon, which lists all the combined notifications
  in one menu. Off by default.
- **Blink**: the tray icon of this kind of element blinks while there are
  unread elements. It only applies to an icon that is not combined. Off by
  default.

Whatever these boxes say, the button of the page in the main window still
shows the count of unread elements.

### Systray tab, group Systray message

- **Group chat**: when checked, a message received in the broadcast chat shows
  a balloon on the RetroShare tray icon. On by default.
- **Chat rooms**: when checked, chat room notifications show a balloon on the
  RetroShare tray icon. Off by default.

When several balloons come in a short time they are grouped into one.

## How to

### Stop all pop-ups for a while

1. Open the **Toasters** tab.
2. Press **Disable All Toasters**. Press it again to get the toasters back.

### Move the toasters away from the task bar

1. Open the **Toasters** tab.
2. Choose a corner in **Position** and adjust **X Margin** and **Y Margin**.
3. Click **Test** to see the result.

## Related settings

- [Sound](settings-sound.md): the sounds played for the same events.
- [Appearance](settings-appearance.md): the toaster button and the sound
  button of the status bar.
- [Chats preferences](settings-chat.md)

## See also

- [Activity](activity.md)
- [Main window, tray icon and status bar](mainwindow.md)
