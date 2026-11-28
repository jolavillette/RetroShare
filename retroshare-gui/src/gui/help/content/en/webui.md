# Web interface

The web interface lets you use your RetroShare node from a web browser. The
node itself serves the pages: you open an address in the browser, log in, and
find your friends, chats, mail, files, channels, forums and boards.

It is useful in two cases:

- next to the desktop program, to reach your node from another device of your
  home network, such as a phone or a tablet;
- with `retroshare-service`, the version of RetroShare that has no window and
  runs on a server or a small computer. The web interface is then the way to
  use the node.

The web interface talks to the node through the JSON API, the programming
interface of RetroShare. It shows the same node as the desktop program: what
you do in one is visible in the other. It does not offer everything the
desktop program does.

Not every RetroShare package contains the web interface. If the tray menu has
no **Show web interface** entry and the preferences have no **Webinterface**
page, your package was built without it.

## Turning it on

1. Open [Preferences > Webinterface](settings-webui.md).
2. Type a password in **Password:**. It protects the access to the web
   interface.
3. Check **Enable Retroshare WEB Interface** and click **Apply settings**.
   The **Status:** light turns green when the web interface is running.

The address and the port on which the node listens are those of the JSON API,
shown in [Preferences > JSON API](settings-jsonapi.md): by default the address
`127.0.0.1`, which only this computer can reach, and the port 9092.

These settings concern the desktop program only. `retroshare-service` has its
own command line options to turn the web interface on.

## Opening it

- Right click the RetroShare tray icon and choose **Show web interface**, or
  click **Start web browser** in the
  [Webinterface preferences](settings-webui.md). Your browser opens the page.
  If the web interface is not turned on, RetroShare tells you so
  (**Webinterface not enabled**).
- Or type the address yourself in a browser: with the default settings,
  `http://127.0.0.1:9092/index.html`.

## Logging in

The login page asks for the password you set in the preferences. The user
name is `webui`; it is already filled in. A link under the password shows more
options: the user name and the address (Url) of the node to connect to, for
the cases where the page is not served by the node itself.

**Log out**, at the bottom of the navigation bar, closes the session.

## Layout

- On a computer screen, a navigation bar on the left lists the sections. It
  can be collapsed to icons. A number next to a section is its count of
  unread items. At the bottom are **Reload**, which reloads the web interface,
  and **Log out**.
- On a phone, the main sections are at the bottom of the screen and the
  others behind a **More** entry.
- A status area shows the number of friends connected out of the total, and
  the state of the connection.

## Sections

### Home

Shows your RetroShare ID with **Copy RetroShare ID** and **Share RetroShare
ID**, and the **Add friend** area: paste the ID you received from a friend,
or drop a text file that contains it, check the details shown and confirm.
The same exchange of IDs is needed as in the desktop program, see
[Adding friends](start-friends.md).

### Network

Lists your friend nodes with their state. Select a friend to see its details
and locations, to chat directly with it (text, emoji, pictures, links to
files) and to remove it (**Remove Friend**). You also set your own status
(**Online**, **Away**, **Busy**) and your status message here. A graph view
draws your network of friends.

### People

Lists identities: your contacts (**Contacts**), your own identities (**My
Identities**) and everyone your node knows (**All Users**). For an identity
you can see its details, **Add to Contacts** or **Remove from Contacts**,
**Start chat** (a distant chat, through an anonymous tunnel) and **Send
mail**.

### Chat

The chat rooms. You can search the rooms, join one, read and write in it, see
the **Participants**, invite friends (**Invite friends to this room**) and
leave it (**Leave Room**). **Create New Chat Room** creates a room with a
name and a topic. Invitations to rooms and unread messages are counted next to
the section name.

### Mail

Your mail, in the folders **Inbox**, **Sent**, **Drafts**, **Outbox** and
**Trash**, with the tags **Important**, **Work**, **Personal**, **Todo** and
**Later** and a filter by tag. Select a message to read it; **Reply**, **Reply
All**, **Forward**, **Delete** and **Mark as unread** act on it. **Compose
Mail** writes a new message, with attachments.

### Files

- **files**: the transfers in progress, downloads and uploads, with their
  size, speed, sources and remaining time. You can add a download from a file
  link (**Add new file**, **Enter the file link:**), cancel a download and
  **Clear completed**.
- **search**: searches files by keywords through the anonymous tunnels, and
  downloads the results you choose.
- **MyFiles**: the directories you share. **Add Shared Directory** shares a
  new one (you type its path on the computer that runs the node) and **Manage
  Visibility** sets who may browse it.
- **FriendsFiles**: the directories your friends let you browse.

### Channels

The tabs **MyChannels**, **Subscribed**, **Popular** and **Other** list the
channels you own, those you follow, the most popular ones and the rest. Open
a channel to subscribe to it, to read its posts, download their attached
files, vote and comment. **Create Channel** creates a channel; in one of your
channels, you can publish a post.

### Forums

The same four tabs for forums: **MyForums**, **Subscribed**, **Popular** and
**Other**. Open a forum to subscribe, read the threads and write messages and
replies. **Create Forum** creates a forum.

### Boards

The same four tabs for boards: **MyBoards**, **Subscribed**, **Popular** and
**Other**. Open a board to subscribe, read the posts page by page, in card or
compact view, vote with the identity you choose and comment. **Create a
Post** adds a post and **Create Board** creates a board.

### Statistics

**Traffic** shows how the data received and sent since the node started is
split between your friends and between the services. **Bandwidth** shows the
current rates and queues per friend, like the Bandwidth page of the
[Statistics](statistics.md) window.

### Config

The settings of the node, in tabs: **network** (addresses, network mode,
dynamic DNS, hidden service), **node** (name, version and public information
of your node), **services** (the services and which are allowed by default),
**files** (download and partial files directories, shared directories,
transfer options), **people** (reputation), **chat** (history, default
identity for chat rooms, who may start a distant chat) and **mail** (tags).
They are the settings of the node itself: a change made here is also seen by
the desktop program.

### Debug

Technical information for diagnosis: the versions of the web interface and of
the node, the response time of the node, the state of the connection between
the page and the node.

## Security

The web interface gives whoever logs in the control of your node: your
messages, your files, your friends. Keep these points in mind:

- **The connection is not encrypted.** The pages and the password travel in
  clear between the browser and the node.
- **The listen address decides who can connect.** With the default
  `127.0.0.1`, only a browser running on the same computer can reach the
  web interface. The option **Allow access from all IP addresses (Default:
  localhost only)** of the [Webinterface preferences](settings-webui.md)
  opens it to the other computers of your network.
- **Do not expose it to the Internet.** To reach your node from outside your
  home, go through an encrypted channel that you control, such as an SSH
  tunnel or a VPN, and keep the listen address on `127.0.0.1`.
- **Choose a real password**, different from the one of your profile. Change
  it in the preferences and click **Apply settings**.
- Turn the web interface off when you do not use it.

## Related settings

- [Webinterface](settings-webui.md): enabling, password, access from other
  computers, directory of the web interface files.
- [JSON API](settings-jsonapi.md): listen address, port and access tokens of
  the JSON API.

## See also

- [Main window, tray icon and status bar](mainwindow.md)
- [Command line options](commandline.md)
- [Adding friends](start-friends.md)
