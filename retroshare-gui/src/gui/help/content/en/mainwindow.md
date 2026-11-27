# Main window, tray icon and status bar

The main window holds all the pages of RetroShare. Around the pages you find
the buttons that switch from one page to another, the status bar at the bottom
and, outside the window, the RetroShare icon in the system tray of your
desktop. This page describes these three elements and the pop-up
notifications.

## Layout

- The title bar reads "RetroShare", the version, "a secure decentralized
  communication platform", then the name of your profile and, in brackets, the
  name of your node.
- The page buttons, as a toolbar or as a list on the left side.
- The current page, in the centre.
- The status bar, at the bottom.

## Page buttons

There is one button per page: **Home**, **Network**, **People**, **Chats**,
**Mail**, **Files**, **Channels**, **Forums**, **Boards**, the pages added by
[plugins](plugins.md) such as **FeedReader**, **Activity** and
**Preferences**. Click a button to show its page. RetroShare reopens on the
page that was displayed when you quit.

After the pages comes one action:

- **Quit**: closes RetroShare, after the confirmation described below.

When a page has something new for you (unread mail, new chat messages, new
posts, new entries in Activity...), its button shows the number in brackets
after the name, for example **Mail (2)**, in bold and red. Hover the button to
read what the number counts.

The buttons are shown either as a toolbar or as a list at the left of the
window. Choose with **Main page items:** (**Buttons** or **Item list**) in the
[Appearance preferences](settings-appearance.md), where you also set the
style of the buttons (icon, text or both) and the size of the icons. In list
mode, **Quit** is the last item of the list, and you can drag the border
between the list and the page to change the width of the list.

### Keyboard

- F1 opens this manual on the topic of the page that is displayed.

## Closing and quitting

- The close button of the window hides the window when **Minimize to Tray
  Icon** is checked in the [General preferences](settings-general.md).
  RetroShare keeps running and stays connected; bring the window back from the
  tray icon. When the option is not checked, the close button quits
  RetroShare.
- **Quit**, in the page buttons or in the tray menu, asks **Do you really want
  to exit RetroShare ?** in a dialog titled **Really quit ?**. Answer **Yes**
  to quit. Check **Do not show the Quit RetroShare MessageBox** in the
  [General preferences](settings-general.md) to quit without the question.

## Tray icon

RetroShare puts its icon in the system tray (the notification area of your
desktop).

### What the icon shows

- The icon has four drawings, according to the number of friend nodes
  currently connected: none, one, two, three or more.
- A star is added on the icon when something new waits for you and the
  notifications are combined in the tray menu (see the **Notify** submenu
  below).
- The tooltip shows "RetroShare", the current speeds (**Down** and **Up**, in
  kB/s), the number of friends connected, the name of your profile and node
  and, when there is something new, the number of new messages. With
  **Disable SysTray ToolTip** checked in the
  [Appearance preferences](settings-appearance.md) the tooltip is not
  refreshed any more.

### Clicking the icon

- On Windows and Linux, a click or a double click shows the main window when
  it is hidden or minimised, and hides it when it is visible.
- On macOS, a double click does the same. A click on the Dock icon shows the
  window.
- A right click opens the tray menu.

### Tray menu

- **Show** or **Hide**: shows or hides the main window. The label depends on
  whether the window is currently visible.
- **Status**: your presence status, as your friends see it. See the next
  section.
- **Notify**: shown only when something new waits for you. It has one entry
  per page concerned, for example "2 new messages". Click an entry to open
  that page.
- **Open Messages**: shows the main window on the [Mail](mail.md) page.
- **Show web interface**: opens the [web interface](webui.md) in your
  browser.
- **Bandwidth Graph**: opens the small window with the curves of your upload
  and download speed. See [Statistics](statistics.md).
- **Statistics**: opens the [Statistics](statistics.md) window.
- **Options**: shows the main window on the [Preferences](settings.md) page.
- **Help**: opens this manual.
- **Minimize**: reduces the main window to the task bar.
- **Maximize**: shows the main window on the whole screen.
- **Quit**: closes RetroShare, after the confirmation described above.

### Status submenu

- **Online**: you are available. This is the normal status.
- **Busy**: tells your friends that you do not want to be disturbed.
- **Away**: tells your friends that you are not at your computer.

The status is only an information shown to your friends next to your name. It
changes nothing in the way your node works: transfers and messages go on.

While your status is **Online**, RetroShare changes it by itself to **Idle**
when you have not touched the keyboard or the mouse for the **Idle Time** set
in the [General preferences](settings-general.md), and back to **Online** when
you return. An **Idle Time** of 0 turns this off.

The same choice is available in the status bar.

### Separate notification icons

In the [Notify preferences](settings-notify.md), tab **Systray**, each kind of
notification can be turned on or off, and can either be combined in the
**Notify** submenu or get its own icon in the system tray. A separate icon is
shown only while something new waits; its tooltip gives the count ("You have 2
new messages") and a click on it opens the page concerned. It can also blink.

### Tray messages

When the options **Group chat** and **Chat rooms** of the group **Systray
message** are checked in the [Notify preferences](settings-notify.md), a new
message in the group chat or in a chat room is announced by a balloon next to
the tray icon. When several messages arrive in a short time, one balloon
announces them together ("3 new messages").

### macOS menus

On macOS, the menu of the Dock icon offers **Open Messages**, **Bandwidth
Graph**, **Statistics**, **Options**, **Help** and the **Status** submenu.
The menu bar adds **Minimize** (Cmd+M) and **Close window**.

## Status bar

The status bar is the line at the bottom of the main window. Each of its
elements can be hidden, and the whole bar too, in the group **Status Bar** of
the [Appearance preferences](settings-appearance.md). **Compact Mode** in the
same group removes the texts and keeps the icons and figures.

From left to right:

- **Status list**: **Online**, **Busy** or **Away**, the same choice as the
  **Status** submenu of the tray menu.
- **Friends**: the number of friend nodes connected, then the total number of
  your friend nodes (tooltip **Online Friends/Total Friends**). The icon is
  grey when nobody is connected.
- **NAT**: a coloured dot that says how reachable your node is from the
  Internet. Hover it to read the state:
    - green, **RetroShare Server** or **Forwarded Port**: other nodes can
      reach you directly;
    - yellow, **Behind Firewall**: you are behind a router, connections are
      still possible with the help of the DHT;
    - yellow, **DHT Disabled**, **Network Restarting** or **Network Status
      Unknown**;
    - red, **Nasty Firewall** or **DHT Disabled and Firewalled**: your friends
      will have trouble reaching you; open a port in your router, see the
      [Network preferences](settings-network.md);
    - grey, **Offline**: no network.
- **DHT**: a coloured dot followed by two figures. The DHT is the public
  directory that your friends' nodes use to find your current address. Hover
  the dot to read the state:
    - grey, **DHT Off**: the DHT is turned off;
    - red, **No peer found in DHT**: the DHT does not answer yet;
    - yellow, **DHT Searching for RetroShare Peers**: the DHT works, few
      RetroShare nodes are found so far;
    - green, **DHT Good**.
  The figures are the estimated number of RetroShare users in the DHT and, in
  brackets, the total number of DHT users.
- **Tor** or **I2P** (hidden nodes only, in place of NAT and DHT): the state
  of the proxy through which your node connects, for example **Tor proxy is
  OK**, **Tor proxy is not available**, **i2p proxy is OK**; with the bundled
  Tor, **Tor is OK** or **Tor is currently offline**.
- **Hashing**: while RetroShare computes the fingerprint of your shared
  files, the text **Hashing file** followed by the progress and the file name.
  Click it to pause the hashing; the text becomes **[Hashing is paused]**,
  click again to resume. **Saving file index...** appears while the list of
  your shared files is written to disk.
- **Discovery**: two counters with an up and a down arrow, the outgoing and
  incoming discovery exchanges that are waiting (tooltips **Waiting outgoing
  discovery operations** and **Waiting incoming discovery operations**).
  Discovery is the service by which your friends tell you about their other
  nodes and their friends.
- **Down** and **Up**: your current download and upload speed in kB/s, each
  followed in brackets by the total amount since RetroShare started. The
  arrows of the icon light up when data flows.
- **Operating mode** (hidden by default): a list that changes the behaviour
  of your node in one click, see below.
- **Sound**: a loudspeaker button. Click it to mute or restore all the sounds
  of RetroShare (tooltips **Sound is on, click to turn it off** and **Sound is
  off, click to turn it on**).
- **Toasters**: a button that turns all the pop-up notifications off or on
  again (tooltips **Toasters are enabled** and **All Toasters are disabled**).
- **Tray icon copy** (hidden by default): a button with the same drawing as
  the tray icon. Click it to open the tray menu. Useful on a desktop that has
  no system tray.

A click on the **NAT** or the **DHT** indicator opens the
[Statistics](statistics.md) window. A click on the rates opens the Bandwidth
Graph.

### Operating mode

- **Normal Mode**: everything is on.
- **No Anon D/L**: turns anonymous tunnels off. Your node no longer downloads
  anonymously and no longer relays the tunnels of others; transfers with your
  direct friends still work.
- **Gaming Mode**: limits RetroShare to a quarter of the speed limits you have
  set, to leave the connection to another program.
- **Low Traffic**: limits RetroShare to a tenth of your speed limits and turns
  anonymous tunnels off.

The first time you choose **No Anon D/L** or **Low Traffic**, a warning
titled **Turtle routing disabled!** reminds you that without tunnels you
cannot download files anonymously nor use distant chat, and that mail is
slower. The mode is kept when RetroShare restarts: do not forget that you
changed it.

## Toasters

A toaster is a small pop-up window that appears in a corner of the screen when
something happens, and goes away by itself. Each one has a close button, and a
main button that takes you to the subject:

- **Friend Connect**: a friend node has just connected. The button opens a
  private chat with this friend.
- **New Message**: a mail has arrived; the toaster shows the sender and the
  subject. The button opens the [Mail](mail.md) page.
- **Download completed**: a download is finished. The button opens the file
  with the program your system associates with it.
- **Private Chat**: a friend has written to you. The button opens the chat
  window.
- **Group Chat**: someone has written in the group chat (broadcast). The
  button opens the group chat in the [Network](network.md) page.
- **Chat Room**: someone has written in a chat room you have joined. The
  button opens the room.
- **Friend Requests**: a node that is not your friend tries to connect to you,
  usually because its owner has added your ID. The button opens the friend
  wizard so that you can accept the request.

Choose which toasters you want, their corner of the screen and their margins
in the [Notify preferences](settings-notify.md), tab **Toasters**. The
**Toasters** button of the status bar and the **Disable All Toasters** button
of that tab silence them all for the moment without changing your choices.

## How to

### Keep RetroShare running in the background

1. Open [Preferences > General](settings-general.md).
2. Check **Minimize to Tray Icon**.
3. Close the main window: RetroShare stays in the system tray.
4. Click the tray icon to bring the window back.

### Stop the notifications for a while

1. In the status bar, click the **Toasters** button to turn the pop-ups off.
2. Click the **Sound** button to mute the sounds.
3. Click both again when you want them back.

## Related settings

- [Appearance](settings-appearance.md): toolbar or list, button style, icon
  size, content of the status bar.
- [Notify](settings-notify.md): toasters, notification icons, tray messages.
- [General](settings-general.md): closing to the tray, quit confirmation,
  idle time, starting minimised.
- [Sound](settings-sound.md): the sound played for each event.

## See also

- [Statistics](statistics.md)
- [Command line options](commandline.md)
- [RetroShare links](links.md)
