# Network

The Network page shows your friends: the profiles you accept connections
with, and their nodes. From here you look at the state of each connection,
arrange your friends in groups, remove a friend, and reach the details of a
node or of a profile.

The page also holds the broadcast chat, a graph of the network around you and
the keyring, the list of all the profile keys your node knows.

## Layout

- At the top left, your avatar, your profile name followed by your node name
  in parentheses, and your status message.
- Below, the search field and the friend list.
- On the right, three tabs: **Broadcast**, **Network graph** and **Keyring**.

The bar between the two sides can be dragged to resize them.

### The friend list

The list is a tree with up to three levels:

- **Groups**, when they are displayed. The name of a group is followed by the
  number of friends connected and the number of friends in the group. A
  friend can be in several groups. Friends that are in no group are listed
  at the top level.
- **Friends**, one row per profile, with the avatar of one of the nodes.
- **Nodes**, one row per node of the profile. A friend who runs RetroShare on
  two computers has two nodes.

Friends and nodes that are connected are written in bold, and the colour of
the name follows the status.

The columns are **Name**, **Last contact** (the last time the node was
connected), **IP** (the address of the connection, `---` when the node is not
connected) and **Id** (the node ID, the profile ID or the group ID). Only
**Name** is displayed at first. Click a column header to sort the list.
Groups always come first.

### Broadcast

The broadcast chat sends your message to all your connected friends at once,
and shows what your friends broadcast. Only your direct friends receive it.
The chat area works like the other chats, see [Chats](chats.md).

### Network graph

The graph draws the profiles around you as discs linked by lines: a line
between two profiles means that they are friends, as far as your node has
learnt it through discovery. Your own profile is yellow, your friends are
green, the friends of your friends are cyan.

- Drag a disc to move it. While you hold a disc, it turns yellow, its
  neighbours green and the others cyan, which shows who is friend with whom.
- Hover a disc to see the name and the profile ID.
- Turn the mouse wheel, or press `+` and `-`, to zoom. Press the space bar
  to shuffle the discs.

### Keyring

The keyring lists every profile key your node knows: yours, those of your
friends, and those received from your friends through discovery. The columns
are:

- **Connections**: **Accepted** when you accept connections with this
  profile, **Denied** otherwise.
- **Profile**: the name of the profile.
- **Trust level**: **Personal signature** when you have signed the key,
  otherwise **Fully trusted peer**, **Marginally trusted peer** or
  **Untrusted peer**, from the trust level you have set.
- **Has signed your key?**: **Yes** or **No**.
- **Id**: the ID of the profile key.
- **Last used**: when the key was used last, to check a signature or a
  connection: **Last hour**, **Today**, a number of days, or **Never**.

Double-click a row to open the [profile details](network-details.md).
Several rows can be selected.

## Buttons

- **?** (help button, top right): opens the manual on this page. The first
  time you open the page, a short introduction to the page is shown in a
  bubble; click the bubble to close it.
- Your avatar: click it to open the **Change Avatar** dialog.
- Your status message (tooltip **Edit status message**): click it to open the
  **Status message** dialog. It shows **Set your status message here.** when
  the message is empty.
- **Search** field above the friend list: keeps the rows that contain all
  the words typed, in the name, and in the IP or the ID when these columns
  are displayed. The tree is expanded while a filter is set, and collapsed
  when you clear it.
- The Delete key removes the selected row, like **Remove Group**, **Deny
  connections** or **Remove Friend Node** below.

In the **Network graph** tab:

- **Redraw**: forgets the positions of the discs and draws the graph again.
- **Friendship level:** how far the graph goes: 0 for yourself only, 1 for
  your friends, 2 to add the friends of your friends.
- **Edge length:** the length of the lines, to spread or tighten the graph.
- The text field: type a part of a name or of a profile ID. The discs that
  match are drawn bigger.
- **Freeze**: stops the automatic layout of the discs.

In the **Keyring** tab:

- **Filter:** keeps the rows that contain the text in any column.
- **Trusted keys only**: keeps the keys whose validity is at least marginal.

## Context menus

### Friend list, on a group

- **Edit Group**: opens the **Edit Group** dialog to rename the group and
  choose its members.
- **Remove Group**: deletes the group at once, without confirmation. The
  friends in it are not removed. Disabled for the standard groups
  (**Friends**, **Family**, **Co-Workers**, **Other Contacts**,
  **Favorites**).

The common entries described below follow.

### Friend list, on a friend

- **Profile details**: opens the [profile details](network-details.md):
  trust level, key signature, options that apply to all the nodes of the
  friend.
- **Deny connections**: after confirmation, stops accepting connections with
  this profile. All its nodes leave the friend list. The key stays in the
  keyring, so you can accept the profile again from the **Keyring** tab.
- **Groups > Create new group**: opens the **Create a Group** dialog. The
  **Groups** submenu is shown only when groups are displayed.
- **Groups > Add to group**: puts the friend in one more group. Shown only
  when you click the friend inside a group, and lists the groups the friend
  is not in.
- **Groups > Move to group**: removes the friend from all groups and puts it
  in the chosen one. Lists the groups the friend is not in.
- **Groups > Remove from group** followed by the name of the group: takes
  the friend out of the group you clicked in. Shown only inside a group.
- **Groups > Remove from all groups**: shown when the friend is in at least
  one group.

The common entries follow.

### Friend list, on a node

- **Chat**: opens a private chat with this node. See [Chats](chats.md).
- **Send message to this node**: opens the mail composer with this node as
  recipient. See [Mail](mail.md).
- **Node details**: opens the [node details](network-details.md): state,
  addresses, certificate.
- **Recommend this node to...**: opens the mail composer with a message that
  contains the certificate link of this node. Choose the friends who receive
  it.
- **Attempt to connect**: starts a connection attempt now and opens the
  **Connection Progress** window, see
  [Adding friends](start-friends.md). Not shown when your node is hidden and
  the other one is not.
- **Copy certificate link**: copies a `retroshare://certificate` link of this
  node to the clipboard. See [RetroShare links](links.md).
- **Remove Friend Node**: after confirmation, removes this node only. The
  other nodes of the friend stay.

The common entries follow.

### Friend list, common entries

These entries end every menu of the friend list. They are the only ones
when you click the empty space with no row selected.

- **Paste certificate link**: opens the friend wizard on the certificate
  link held in the clipboard. Disabled when the clipboard holds no such link.
- **Expand all**: opens all groups and friends.
- **Collapse all**: closes them.
- **export friendlist**: saves your friends, with the certificates of their
  nodes and your groups, to an XML file. The file is not encrypted.
- **import friendlist**: adds the friends and groups of such a file. A
  message tells whether some friends or group memberships could not be added.

### Friend list, on the column headers

- **Online friends on top**: lists the connected friends before the others,
  whatever the sort column.
- **Offline Friends**: shows the friends that are not connected. Unchecked,
  only the connected ones are listed.
- **Status**: writes a second line under each node name, with the status
  message of the node, or its status in parentheses when it has no message.
- **Status icons**: adds the status icon of the friend on its avatar.
- **Groups**: shows the groups. Unchecked, the list holds the friends only.
- **Last contact**, **IP**, **ID**: show or hide these columns.

### Network graph, on a disc

- **Deny friend** (on a friend) or **Make friend** (on another profile): both
  open the [profile details](network-details.md), where the **Deny
  connections** or **Accept connections** button does it. Not shown on your
  own disc.
- **Peer details**: opens the profile details.

### Keyring, on a key

- **Profile details...**: opens the [profile details](network-details.md) of
  the key. From there you can accept connections with a profile that is not
  your friend.
- **Remove unused keys...**: opens the **Clean keyring** dialog. It lists the
  keys, with those that have not been used for three months and are not
  accepted already checked. The keys checked are deleted when you confirm.
- **Remove this key**: deletes the selected keys from the keyring. If some of
  them are friends, a confirmation is asked and they are removed from your
  friends first.

Before any deletion the keyring is saved to a backup file, whose name is
shown. Your own key cannot be deleted.

## Dialogs

### Create a Group, Edit Group

- **Group Name:** the name of the group. **OK** is disabled while the name
  is empty or already used. The name of a standard group cannot be changed.
- **Group ID:** the identifier of the group, **To be defined** until the
  group is created.
- **Friends**: check the friends that are members of the group. The search
  field filters the list.

Groups are also used to restrict who sees your shared files, see
[Shared files](files-shared.md).

### Change Avatar

- **Browse...**: loads a picture file (JPG or PNG). Use the mouse to zoom
  and move the picture in the frame.
- **Remove**: removes your avatar.
- The stickers under the picture, when you have some, can be used as avatar
  with one click.

The avatar is the one of your node: your friends see it in their friend
list and in private chats. It is not the avatar of your identities.

### Status message

- **Message:** a free text that your friends see next to your node, for
  example in their friend list when **Status** is checked.

## How to

### See who is connected

1. Right-click the column header of the friend list.
2. Uncheck **Offline Friends**, or check **Online friends on top**.

### Put friends in a group

1. Right-click a friend and choose **Groups > Create new group**.
2. Type a name, check the members and click **OK**.
3. To add someone later, right-click the group and choose **Edit Group**, or
   right-click the friend and use **Groups > Move to group**.

### Remove a friend

1. Right-click the friend (the profile row, not a node).
2. Choose **Deny connections** and confirm.

### Become friend with a profile of the keyring

1. Open the **Keyring** tab and double-click the profile.
2. Click **Accept connections**.
3. Your node needs an address to reach the other one. If the connection does
   not come, ask the person for their RetroShare ID and add it with the
   friend wizard, which is the recommended way.

## Related settings

- [Network](settings-network.md): your addresses, discovery, DHT.
- [Permissions](settings-permissions.md): which services each friend may
  use.
- [Notify](settings-notify.md): notifications when a friend connects.

## See also

- [Friend details](network-details.md)
- [Adding friends](start-friends.md)
- [Home](home.md)
- [People](people.md)
- [Activity](activity.md)
