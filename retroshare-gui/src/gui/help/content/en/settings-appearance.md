# Appearance

This page of the preferences sets how RetroShare looks: the language, the
format of dates, the navigation of the main window, the size of the text, the
contents of the status bar, and the style of the windows. The elements it
changes are described in
[Main window, tray icon and status bar](mainwindow.md).

Every change is applied at once, except the language, which needs a restart.

## Layout

The page is a set of groups: **Language**, **Date Format**, **Tool Bar**,
**Fonts**, **Status Bar**, **Style** and **Style Sheet**.

## Buttons

### Language

- The list (tooltip **Choose the language used in RetroShare**): the language
  of the interface, shown with its flag. The first time, RetroShare takes the
  language of your system. The change is saved at once but the text
  **(Needs restart)** next to the list reminds you that the windows change
  language only after you restart RetroShare.

### Date Format

The list (tooltip **Choose the date display format**) sets how dates and times
are written in the lists of RetroShare. Each choice shows the current date as
an example:

- **System Default**: the short date and time format of your system. This is
  the default.
- **ISO 8601**: year-month-day followed by the time on 24 hours, for instance
  `2026-11-27 14:05`.
- **Text**: the long date format of your system, with the name of the month,
  followed by the time.

### Tool Bar

- **Main page items:** chooses how you go from one page of the main window to
  another:
  - **Buttons**: a tool bar with one button per page. This is the default.
  - **Item list**: a list with one line per page; the tool bar is hidden.
- The list under **Choose the style of Tool Buttons.** sets what a button of
  the tool bar shows: **Icon Only**, **Text Only**, **Text Beside Icon** or
  **Text Under Icon** (the default). It only concerns the tool bar, and is
  hidden when you switch to **Item list**.
- The icon size list sets the size of the icons of the tool bar and of the item
  list: **Icon Size = 8x8**, **16x16**, **24x24**, **32x32**, **64x64** or
  **128x128**. The default is chosen from the size of the font of your system,
  32x32 on an ordinary screen.

### Fonts

- **Minimum font size**: the size, in points, of the text in the main lists of
  RetroShare: the friend list, the lists of people, mails, chat rooms, search
  results, forum threads, the lists of channels, forums and boards, and the
  list of pages of the preferences. From 5 to 64; 11 by default, 13 on macOS.

### Status Bar

The group has a checkbox in its title: uncheck **Status Bar** to hide the whole
status bar at the bottom of the main window. The boxes inside choose its
elements:

- **Show Status ComboBox**: the list where you set your own status. On by
  default.
- **Show Peer Status**: the counter **Friends:** with the number of connected
  friend nodes over the total. On by default.
- **Show NAT Status**: the **NAT:** light, which tells whether your node can be
  reached from the Internet. On by default. Not shown on a hidden node.
- **Show DHT Status**: the **DHT:** light and the number of RetroShare nodes
  found in the DHT. On by default. Not shown on a hidden node.
- **Show Hashing Status**: the text that names the file being hashed, where you
  click to pause or resume hashing. On by default.
- **Show Discovery Status**: the icons that appear while discovery exchanges
  are waiting. On by default.
- **Show Network Rate Status**: the **Down:** and **Up:** transfer rates. On by
  default.
- **Show Operating Mode Status**: the list that switches the node between
  **Normal Mode**, **No Anon D/L**, **Gaming Mode** and **Low Traffic**. Off by
  default.
- **Show Sound Status**: the button that mutes and unmutes all sounds. On by
  default.
- **Show Toaster Disable**: the button that disables and enables all
  toasters. On by default.
- **Show SysTray on Status Bar**: a button with the RetroShare icon that opens
  the menu of the tray icon. Useful on a desktop that has no system tray. Off
  by default.
- **Compact Mode** (tooltip **Remove surplus text in status bar.**): shortens
  the friends, NAT, DHT, hashing and rate elements so that the status bar takes
  less room. Off by default.
- **Disable SysTray ToolTip**: when checked, the tooltip of the tray icon only
  shows "RetroShare" and the name of your node; it no longer shows the
  transfer rates and the number of connected friends. Off by default.

### Style

- The list (tooltip **Choose RetroShare's interface style**): the widget style
  used to draw buttons, lists and windows. The choices are the styles your Qt
  installation provides, for instance Fusion or Windows; they differ from one
  system to another.

### Style Sheet

- The list: the style sheet, which sets the colours and decorations on top of
  the style. RetroShare comes with **Standard_Light (built-in)**, the default,
  **Standard_Dark (built-in)**, the dark theme, and **default (built-in)**. The
  `.qss` files found in the `qss` directory of your RetroShare configuration
  directory or of the RetroShare data directory are added to the list under
  their file name.

## How to

### Switch to the dark theme

1. In **Style Sheet**, choose **Standard_Dark (built-in)**.

### Get more room for the pages

1. In **Tool Bar**, choose **Icon Only** and a smaller icon size, or switch
   **Main page items:** to **Item list**.
2. In **Status Bar**, check **Compact Mode** or uncheck the elements you do not
   use.

## Related settings

- [General](settings-general.md)
- [Notify](settings-notify.md): toasters and tray icons.
- [Sound](settings-sound.md)
- [Forums preferences](settings-forums.md) and
  [Chats preferences](settings-chat.md): font size and contrast of the
  messages.

## See also

- [Main window, tray icon and status bar](mainwindow.md)
