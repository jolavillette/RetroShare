# Activity

The Activity page is a log of what happened on your node since you started
RetroShare: friends that connected, connection attempts by nodes you do not
know, new mail, new channels, forums and boards, new posts and comments, and
circle invitations. Each event is shown as an item, in the order in which your
node received it, and most items carry buttons to react at once.

The page is not a store. It holds the 500 most recent items; older ones are
dropped, and the list is empty again after a restart. Removing an item from
the page never deletes the mail, the post or the friend it refers to.

Which kinds of events are shown is chosen in the
[Notify preferences](settings-notify.md). The number shown on the Activity
button of the main window counts the items of the page, except the friend
connection items.

## Layout

- At the top, the title **Activity Stream**, the sort selector, the **Remove
  All** button and the help button.
- Below, the list of items. Each item has a title line, a row of buttons and,
  for most kinds, a details part that is folded until you expand it.

## Buttons

### Page

- Sort selector: **Newest on top** or **Oldest on top**. The choice is
  remembered.
- **Remove All**: removes every item from the page.
- **?** (help button): opens the manual on this page. The first time you open
  the page, a short introduction to the page is shown in a bubble; click the
  bubble to close it.

### Buttons found on every item

- Arrow button (tooltips **Expand** and **Hide**): unfolds or folds the
  details of the item.
- Cross button (tooltip **Remove Item**): removes this item from the page.
- Tick button (tooltip **Set as read and remove item**), on posts and
  comments: marks the post or comment as read, then removes the item.

Names of channels, forums, boards, posts and people that appear as links in
an item open the corresponding page when you click them.

### Friend items

Shown when a friend node connects (**Friend Connected**), when a peer refuses
your connection (**Connection refused by peer**) and when the clock of a
friend differs from yours (**Friend Time Offset**).

- **Start Chat** (tooltip **Chat**): opens a private chat with the friend
  node. Enabled while the friend is connected.
- **Send Message**: opens the mail composer with the friend node as
  recipient. Enabled while the friend is connected.

The details give **Name:**, **Peer ID:**, **Trust:**, **Location:**, **IP
Address:**, **Connection Method:**, **Status:** and **Time offset:**.

### Security items

Shown when a connection is attempted and something is wrong or unknown:

- **Connect Attempt** and **Peer/node not in friendlist**: a RetroShare node
  that is not your friend tried to connect. This usually means that someone
  added your ID and waits for you to add theirs.
- **Certificate has wrong signature!! This peer is not who he claims to be.**
- **Connection attempt from ... Not a Retroshare user.**: the node presented
  no usable certificate.
- **Unknown Security Issue**: any other case, for example a connection from
  a blacklisted IP address.

The buttons depend on what is known about the node:

- **Accept Friend Request**: shown for a node that is not your friend. Opens
  the friend wizard on this node, so that you can add it as a friend. See
  [Adding friends](start-friends.md).
- **Peer details**: opens the details of the node, or of its profile when the
  node itself is unknown. See [Friend details](network-details.md).
- **Start Chat**: shown when the node is a connected friend.
- **Deny friend**: shown when the node belongs to one of your friends.
  Removes this friend after confirmation.

With the matching option of the [Notify preferences](settings-notify.md), a
connection attempt also puts a notification message in your Inbox.

### Wrong IP address items

**Wrong external ip address reported**: a friend reports seeing you under
another address than the one your node believes it has. This is harmless
just after your address changed. Otherwise it can mean that your connection
to this friend is being relayed by a third party.

- **Peer details**: opens the details of the friend node.
- **Add IP to whitelist**: shown when the address can be whitelisted. Its menu
  lets you add **Only IP** or an **Entire range**; a line then confirms that
  the address was added.

### Non-RetroShare connection items

**Non-Retroshare connection attempt (TLS)**: something that is not RetroShare
connected to your port, most of the time a bot scanning the Internet. These
attempts are collected in a single item; expand it to see the list of **Origin
IP Address** and **Time**. There is nothing to do about them.

### Mail items

**Message From**: a new mail arrived. The item shows the sender and the
subject; expand it to read the text and see the attached files.

- **Reply Message** (tooltip **Reply to Message**): opens the composer on a
  reply.
- Trash button (tooltip **Delete Message**): deletes the mail itself, for
  good, and removes the item.
- **Play Media**: has no effect in this version.
- **Send Invite**: only on a friend request sent by an identity. Opens the
  composer on an answer that contains a link to your own certificate. On such
  a request the reply button reads **Reply to invite**.

Each attached file has its own line with **Download** (tooltip **Download
File**), a cancel button (tooltip **Cancel Download**), **Play** (tooltip
**Play File**, enabled once you have the file) and a button with the tooltip
**Copy RetroShare Link**.

### Channel items

- **New Channel:**: a channel that you can now subscribe to, or whose
  description changed. **Subscribe** (tooltip **Subscribe this Channel**)
  subscribes you; it is disabled when you are already subscribed. The link
  button (tooltip **Copy RetroShare Link**) copies a link to the channel.
  The details show the description.
- **Publish permission received for channel:**: the owner of a channel gave
  you the right to post in it.
- **Channel Feed:**: a new post in a channel. The item shows the channel, the
  title of the post and, for a post with files, their number and total size.
  **Download** fetches all the files of the post. **Play**, or **Open** when
  the files are not media files, opens the files that are complete. Both
  buttons are shown only when the post has files. Expand the item to read the
  whole post and to see each file with its own buttons.
- **New Comment**: a comment on a channel post, with its author and the title
  of the post. The link button (tooltip **Copy RetroShare Link**) copies a
  link to the post.

### Forum items

- **New Forum** and **Forum updated**: a forum you can subscribe to, or whose
  description changed. **Subscribe** (tooltip **Subscribe to Forum**) and the
  link button work as for channels.
- **Moderator list changed**: the details list the moderators that were added
  and removed.
- **Forum Feed:**: a new message in a forum, with its author and subject.
  Expand it to read the message and, for an answer, the message it replies
  to. Click the forum name or the subject to open it in [Forums](forums.md);
  this also removes the item.

### Board items

- **New Board**: a board you can subscribe to. **Subscribe** (tooltip
  **Subscribe to Posted**) and the link button work as for channels.
- **Board:** followed by **Post:**: a new post on a board. Click its picture
  to see it full size.
- **New Comment**: a comment on a board post. The link button copies a link to
  the post.

### Circle items

- **You received a membership request a circle you're administrating:**:
  **Grant membership** adds the identity to the circle.
- **You received an invitation to join this circle:**: **Accept** makes your
  identity a member.
- "... which you invited, has joined this circle you're administrating.":
  **Revoke membership** removes the identity from the circle.
- "... has joined this circle." and "... has left this circle.": for your
  information.
- **Your identity ... as been accepted in this circle.**: the **Cancel
  membership** button of this item has no effect in this version; leave the
  circle from the [People](people.md) page.
- **Your identity ... has been revoked from this circle.**: **Cancel
  membership request** withdraws your request.

On a circle item the arrow button (tooltip **Details**) opens the details
window of the circle instead of unfolding the item. See [People](people.md).

## How to

### Accept someone who added you

1. Look for a **Connect Attempt** item.
2. Click **Peer details** and check that the profile is the one you expect.
3. Click **Accept Friend Request** and finish the wizard.

### Catch up after a while

1. Choose **Newest on top**.
2. Read the items; use the tick button on the posts you have seen.
3. Click **Remove All** when you are done.

## Related settings

- [Notify](settings-notify.md): the kinds of events shown on this page, and
  the pop-ups and tray notifications.
- [Sound](settings-sound.md): sounds played on events.

## See also

- [Network](network.md)
- [Mail](mail.md)
- [Channels](channels.md), [Forums](forums.md), [Boards](boards.md)
- [People](people.md): circles.
