# Chats

The Chats page of the [Preferences](settings.md) sets how chat messages are
written, displayed, searched and kept, for chat rooms, private chats and
distant chats.

## Layout

The page has four tabs: **General**, **Chat rooms**, **Private chat** and
**History**.

## Buttons

### General tab, General settings

- **Minimum text contrast**: from 1 to 21, 4.5 by default. When the colour of
  a text is too close to the background, RetroShare changes it until this
  contrast is reached, so that a message stays readable. Raise the value if
  some messages are hard to read with your style.
- **Send message with Ctrl+Return**: the message is sent with Ctrl+Return and
  the Return key alone starts a new line. When the box is not checked, Return
  sends the message. Off by default.
- **Send as plain text by default**: your messages are sent without
  formatting. You can still change it for one chat in the menu of its window.
  Off by default.
- **Load embedded images**: the images contained in the messages are
  displayed. On by default.
- **Shrink text edit field when not needed**: in a chat room, the field where
  you type is reduced when you are reading the messages, and gets its size
  back when you type. Off by default.

### General tab, Search

These options set the search bar of the chat windows. They are read when a
chat is opened.

- **Show Bar by default**: this box is not saved and has no effect in this
  version.
- **Threshold for automatic search**: the number of characters you must type
  before the search starts by itself. 4 by default, 100 at most. Pressing
  Return starts the search whatever the length.
- **Case sensitive**: the search distinguishes capital and small letters.
- **Whole Words**: the search only finds whole words.
- **Move to cursor**: the conversation scrolls to the text that was found. On
  by default.
- **Color All Text Found**: every occurrence is highlighted, without limit.
  This needs more processing on a long conversation.
- **Maximum count for coloring matching text**: the number of occurrences
  that are highlighted when the previous box is not checked. 40 by default,
  1000 at most.
- **Color of found text**: the button shows the highlight colour; click it
  to choose another one.

### General tab, Fonts

- **Change Chat Font**: opens the font chooser to pick the default font of
  chat messages. **Chat Font:** shows the name of the current font, written
  in that font.
- **Enable custom fonts**, **Enable custom font size**, **Enable custom font
  color**, **Enable bold**, **Enable italics**: each box lets the messages you
  receive keep that part of their formatting. When a box is not checked, that
  formatting is removed and the message is shown with your own font, size,
  colour or weight. All are on by default.
- **Minimum displayed font size**: text written smaller than this size is
  enlarged to it. From 1 to 64, 10 by default.

At the bottom of the tab, two squares help you check the size settings of
your system: the first one measures 1 cm when the resolution your system
reports is right, the second one is scaled with the font size of your
system. Hover them for the figures.

### Chat rooms tab

- **Default identity for chat rooms:** the identity that is proposed when
  you enter a chat room.
- **Blink tab icon**: the icon of the tab of a chat room blinks when a new
  message arrives in it.
- **Hide join/leave messages**: the lines saying that someone entered or left
  the chat room are not shown.
- **Count all unread messages**: every unread message of your chat rooms is
  counted in the notification of new chat room messages.
- **Count occurrences of my current identity**: the messages that contain the
  name of the identity you use in the chat room are counted.
- **Count occurrences of any of the following texts (separate by
  newlines):** the messages that contain one of the texts typed in the box
  below are counted. Write one text per line. The box is greyed when the
  option is not checked.
- **Case sensitive search**: the name of your identity and the texts above
  must appear with the same capital and small letters to be counted.
- **Style:**, **Variant:**: the style used to display the messages of chat
  rooms. See "Styles" below.

### Private chat tab

- **Accept chat from:** who may open a distant chat with you, that is a chat
  through anonymous tunnels with one of your identities:
  - **Everyone**: any identity.
  - **Contacts only**: only the identities that are in your contacts.
  - **Nobody**: distant chats are refused.
- **Open Window for new chat**: a chat window opens when a friend starts a
  chat with you. When the box is not checked, you are only notified and you
  open the chat yourself. A distant chat always opens a window.
- **Grab Focus when chat arrives**: the chat window comes to the front and
  takes the keyboard when a message arrives.
- **Use a single tabbed window**: all the private chats share one window,
  with one tab per chat. When the box is not checked, each chat has its own
  window.
- **Blink window/tab icon**: the icon of the chat window or of its tab blinks
  when a new message arrives.
- **Enable Emoticons**: the text smileys of private chat messages are shown
  as pictures. On by default.
- **Style:**, **Variant:**: the style used to display private chats. See
  "Styles" below.

### History tab

The table at the top has one column per kind of chat: **Chat rooms**,
**Node-to-node chat** (the chat with a friend node), **Broadcast** and
**Distant chat**.

- **Enabled:** the messages of this kind of chat are saved on the disk.
- **Saved messages (0 = unlimited):** the number of messages that are kept
  for each conversation of this kind.
- **Number of messages restored (0 = off):** the number of saved messages
  that are shown again when you open a chat of this kind. By default 20 for
  node-to-node and distant chats, 0 for chat rooms and broadcast.

Below the table, in the **General** group:

- **Maximum storage period, in days (0=keep all):** the saved messages older
  than this are deleted. It keeps the history from filling up with short
  lived conversations, such as those of chat rooms and distant chats.

At the bottom, **Style:** and **Variant:** set the style of the history
window.

### Styles

The **Chat rooms**, **Private chat** and **History** tabs each end with the
same block:

- **Style:** the list of the styles installed for this kind of display.
- **Variant:** the variants of the selected style, for instance other
  colours. The list is greyed when the style has no variant.
- **Author:** and **Description:** are given by the style.
- The preview shows sample incoming, outgoing and system messages with the
  selected style and variant.

A new style is applied to the open chats at once.

## How to

### Keep the history of your private chats for one month

1. Open the **History** tab.
2. Check **Enabled:** in the **Node-to-node chat** column.
3. Set **Saved messages (0 = unlimited):** to 0 in the same column.
4. Set **Maximum storage period, in days (0=keep all):** to 30.

## See also

- [Chats](chats.md)
- [People](people.md)
- [Notify](settings-notify.md)
- [Appearance](settings-appearance.md)
