# Sound

RetroShare can play a sound when something happens: a friend connects, a
message arrives, a download finishes. This page of the preferences chooses
which events make a sound and which file is played.

Every change is applied at once.

## Layout

- At the top, a frame shows the event selected in the list: its name after
  **Event:** and the sound file after **Filename:**, with the buttons
  **Default**, **Browse** and the play button. The frame is greyed until you
  select an event.
- Below, the list of events has two columns: **Event**, with a checkbox, and
  **Filename**. The events are grouped:
  - **Friend** > **Go Online**: a friend node connects.
  - **Chatmessage** > **New Msg**: a chat message arrives.
  - **Message** > **Message arrived**: a new mail arrives.
  - **Download** > **Download complete**: a download finishes.
  - **Lobby** > **Message arrived**: a message arrives in a chat room.
  - A plugin can add its own group of events, for instance the ring of an
    incoming call.

## Buttons

- The checkbox of an event: when checked, the sound of this event is played.
  When unchecked, the event is silent. An event that is checked but has no file
  is silent too.
- **Filename:** field: the path of the sound file of the selected event. You
  can type it or paste it; the cross in the field clears it.
- **Default**: puts back the sound file that comes with RetroShare for this
  event. The button is greyed when RetroShare has no default file for the
  event, or when that file is not installed.
- **Browse**: opens a file dialog, **Open File**, to choose the sound file.
  Only `.wav` files are offered.
- Play button (the arrow): plays the file of the selected event at once, so
  that you can hear it. It is greyed when the event has no file. It plays even
  when the event is unchecked or the sound is muted.

## Default sounds

The first time RetroShare runs, every event whose default sound file is
installed is enabled with that file. The default files are in the `sounds`
directory of the RetroShare data directory.

## Muting all sounds

The sound button of the status bar switches all the sounds off and on without
changing this page. Its tooltip reads **Sound is on, click to turn it off** or
**Sound is off, click to turn it on**. See
[Main window, tray icon and status bar](mainwindow.md).

## How to

### Change the sound of an event

1. Select the event in the list.
2. Click **Browse** and choose a `.wav` file.
3. Click the play button to check it.
4. Make sure the box of the event is checked.

### Silence one event

1. Uncheck the box of the event in the list.

## Related settings

- [Notify](settings-notify.md): the visual notifications of the same events.
- [Appearance](settings-appearance.md): **Show Sound Status** shows or hides
  the sound button of the status bar.
