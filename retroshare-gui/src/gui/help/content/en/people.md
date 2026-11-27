# People

The People page lists the identities your node knows and the circles they
form. An identity is the name under which someone writes in chat rooms,
forums, channels and boards, and the address used for mail and distant chat.
It is not a friend node: you see the identities of people you are not connected
to, and you can talk to them without knowing which node they use.

An identity is either a **pseudonym**, which nobody can link to a node, or
**linked to your profile**, which means it is signed by the profile of its
owner. A signed identity is easier to trust, but it shows which profile wrote
a message. You can own several identities of both kinds.

A **circle** is a group of identities. Circles are used to restrict who can
see and use a forum, a channel or a board.

## Layout

- The title bar shows **People** and, at the right, the **?** help button.
- On the left, a toolbar with the **Search...** field and the **Create new...**
  button, and below it the list of identities.
- On the right, two tabs: **Person**, the details of the identity selected in
  the list, and **Circles**, the list of circles.

### The list of identities

The identities are sorted in three categories:

- **My own identities**: the identities you created.
- **My contacts**: the identities you added to your contacts.
- **All**: every other identity your node has received.

The columns are **Name**, **Id**, **Owner Id**, **Owner** and **Reputation**.
Only the name and the reputation icon are shown at first; right click the
header of the list to show the others. The owner columns are filled for signed
identities only. Hover a line to see who owns the identity. Click a column
header to sort the list.

The **Search...** field keeps the identities whose text matches what you type.
It searches the columns that are visible: the name, and the ID, owner ID and
owner name when you have shown these columns.

### The Person tab

- At the top, the avatar and the name of the identity, with the number of
  **Positive votes** and **Negative votes** given by your friend nodes.
- Under the avatar, the **Send Invite** button, or **Edit Identity** for one of
  your own identities.
- **Identity ID :** the unique identifier of the identity.
- **Owner node ID :** the ID of the profile that signed the identity. It is
  not shown for a pseudonym. `[unverified]` is appended when the signature has
  not been checked yet.
- **Owner node name :** the name of that profile, shown when your node knows
  it.
- **Type:** tells whether the identity is yours, anonymous, or linked to a
  friend node, a known node or an unknown node.
- **Created on :** and **Last used:** the creation date and the last time your
  node saw the identity being used.
- **Your opinion:**, **Friend votes:**, **Overall:** and **Ban-option:** are
  the reputation controls, described below.
- **Usage statistics**: where your node has seen the identity during this
  session, for example as the author of a forum message. The names of forums,
  channels and posts are links that open them.

### Reputation

Reputation decides whether the content signed by an identity is shown,
forwarded or dropped.

- **Your opinion:** is your own vote: **Negative**, **Neutral** or
  **Positive**. A negative opinion bans the identity: its messages are hidden
  and your node stops forwarding them. Your opinion is sent to your friend
  nodes.
- **Friend votes:** counts the positive and negative opinions of your friend
  nodes.
- **Overall:** is the result. Your own opinion wins when it is not neutral:
  **Positive** or **Negative (Banned by you)**. Otherwise the votes of your
  friends decide: **Positive (according to your friends)**, **Negative
  (according to your friends)** or **Neutral**.
- **Auto-Ban all identities signed by the same node**: bans every identity
  signed by the same profile, including those it will create later. The option
  is shown for signed identities only and cannot be changed for your own
  identities.

### How long identities are kept

Your own identities and your contacts are kept for ever. The other identities
are deleted when they have not been used for some time: 5 days for a
pseudonym, 8 days for an identity signed by a profile your node does not know,
30 days for an identity signed by a known profile. Banned identities are
deleted after the delay set in the [People preferences](settings-people.md).
Add someone to your contacts to keep the identity.

### The Circles tab

The circles are listed in two categories, **Circles I belong to** and
**Other circles**. The columns are **Circle name** and **Membership**.

- A circle written in bold is a circle you administrate.
- Open a circle to see the identities attached to it, each with its state in
  the **Membership** column: **Invited**, **Subscription pending** or
  **Member**. Your own identities are written in bold.
- Hover a circle to see its ID, its visibility, your role, and whether you are
  a full member.

## Buttons

- **?** (help button, top right): shows a short introduction to identities and
  circles.
- **Create new...**: opens a menu with **Create new Identity** and **Create
  new circle**. Both open the dialogs described below.
- **Send Invite** (Person tab, not for your own identities): opens a new mail
  to this identity that contains your RetroShare ID and asks the person to
  become your friend. Send it from the mail window. The message stays in your
  outbox until the recipient's node has acknowledged it.
- **Edit Identity** (Person tab, only for your own identities): opens the
  **Edit identity** dialog.
- **Close** (on the yellow note shown after an invite): hides the note.

## Context menus

### An identity in the list

Several identities can be selected; the entries then apply to all of them.

- **Chat with this person**: opens a distant chat with the identity, an
  encrypted conversation that travels through anonymous tunnels. Only when a
  single identity is selected and it is not one of yours. The entry is greyed
  out when you have no identity of your own.
- **Chat with this person as...**: replaces the previous entry when you own
  several identities. The submenu lists them; choose the one you want to talk
  as.
- **Send message**: opens a new mail addressed to the selected identities. Not
  shown when one of your own identities is in the selection. RetroShare asks
  for a confirmation above 20 recipients.
- **Add to Contacts**: moves the identities to **My contacts**. Shown when none
  of the selected identities is a contact yet.
- **Remove from Contacts**: shown when all the selected identities are
  contacts.
- **Copy identity to clipboard**: copies a RetroShare link that contains the
  identity. Paste it in a mail or a chat to hand the identity to someone. Only
  when a single identity is selected.
- **Set positive opinion**, **Set neutral opinion**, **Set negative opinion**:
  set your opinion for all the selected identities. An entry is hidden when
  one of the selected identities already has that opinion.
- **Edit identity**: opens the **Edit identity** dialog. Only for a single
  identity that you own.
- **Delete identity**: deletes one of your own identities after a
  confirmation. This cannot be undone.

### The header of the list

- **Id**, **Owner Id**, **Owner Name**, **Reputation**: show or hide the
  column.

### A circle

- **Edit Circle**: opens the circle dialog, where you change the name, the
  visibility and the list of invited members. Only for a circle you
  administrate.
- **See details**: opens the same dialog, read only. For the circles you do
  not administrate.
- **Copy retroshare link**: copies a RetroShare link to the circle.
- **Request subscription for identity ...**: asks to join the circle with one
  of your identities. You become a member when the administrator grants the
  request.
- **Accept circle invitation for identity ...**: shown when the administrator
  has invited one of your identities. Accepting makes it a full member.
- **Quit this circle for identity ...**: withdraws a member identity from the
  circle.
- **Cancel subscribe request for identity ...**: withdraws a request that has
  not been granted yet.

Each of the four membership entries is shown for the identities it applies
to. When several of your identities are in the same state, the entry becomes
a submenu that lists them.

### An identity inside a circle

- **Request subscription**, **Accept circle invitation**, **Quit this
  circle**, **Cancel subscribe request**: as above, for this identity. Only on
  one of your own identities.
- **Grant membership**: adds the identity to the invited members. Only when
  you administrate the circle, on an identity that asked to join.
- **Revoke this member**: removes the identity from the invited members. Only
  when you administrate the circle.

## Dialogs

### New identity

- **Type**: **Linked to your profile** creates an identity signed by your
  profile; RetroShare asks for your profile password. **Pseudonym** creates an
  identity that cannot be linked to your node.
- **Nickname**: the name other people see, from 2 to 30 characters. Nicknames
  are not unique; the ID tells two identities with the same name apart.
- **Key ID**, **PGP Hash**, **PGP Id**, **PGP Name**: information fields,
  filled by RetroShare.
- **Choose image...**: picks a picture file as the avatar. Use the mouse on
  the picture to zoom and move it until the frame shows what you want.
  Without a picture, RetroShare draws a default avatar from the ID.
- **Remove**: removes the avatar.
- **Create**: creates the identity. **Cancel** closes the dialog.

### Edit identity

The same dialog, for one of your identities. You can change the nickname and
the avatar; the type cannot be changed. **Update** saves the changes.
RetroShare asks for your profile password when the identity is linked to your
profile.

### Create New Circle

- **Name:** the name of the circle.
- **Contact author:** the identity of yours that is published as the contact
  of the circle.
- **Distribution:** who can see the circle and its list of members.
  - **Public**: everybody.
  - **Private**: only the invited members.
  - **Only visible to members of:** only the members of the circle you choose
    in the list next to the option.
- **Invited Members**: the identities you invite. The number of members is
  shown in the title.
- **Known People**: the identities your node knows. **All**, **Signed** and
  **Signed by friend node** choose which ones are listed, and **Filter**
  searches them by nickname or ID.
- **Add**: invites the identity selected in **Known People**. The context menu
  of that list has the same entry, **Add Member**.
- **Remove**: removes the identity selected in **Invited Members**. The
  context menu of that list has the same entry, **Remove Member**.
- **Create**: creates the circle. A message shows its name and ID.

The same dialog, titled **Edit Circle**, changes a circle you administrate;
its button is **Update**. Titled **Circle Details**, it shows a circle you do
not administrate: nothing can be changed and the only button is **Close**.

## How to

### Create an identity

1. Click **Create new...** and choose **Create new Identity**.
2. Choose **Pseudonym**, or **Linked to your profile** if you want other
   people to see that the identity belongs to your profile.
3. Type a nickname, choose an avatar if you want one, and click **Create**.

### Ban someone

1. Select the identity in the list.
2. Right click it and choose **Set negative opinion**, or set **Your
   opinion:** to **Negative** in the **Person** tab.
3. If the identity is signed and you want to ban every identity of the same
   profile, tick **Auto-Ban all identities signed by the same node**.

### Become a member of a circle

An identity is a member of a circle when two things are true: the
administrator has invited it, and the identity has asked to join. The order
does not matter.

1. Open the **Circles** tab and find the circle under **Other circles**.
2. Right click it and choose **Request subscription for identity ...**. The
   identity appears under the circle as **Subscription pending**.
3. Wait for the administrator to grant the request. The state becomes
   **Member** and the circle moves to **Circles I belong to**.

If the administrator invited you first, the identity is shown as **Invited**:
right click and choose **Accept circle invitation**.

### Restrict a forum to a group of people

1. Click **Create new...** and choose **Create new circle**.
2. Name the circle, choose **Private** if the list of members must stay
   hidden, and add the identities to **Invited Members**.
3. Click **Create**. Each invited person must accept the invitation.
4. Create the forum, channel or board and restrict it to the circle. See
   [Forums](forums.md), [Channels](channels.md) and [Boards](boards.md).

## Related settings

- [People](settings-people.md): default opinion for contacts, automatic
  contacts, how long banned identities are kept.

## See also

- [Chats](chats.md)
- [Mail](mail.md)
- [Glossary](glossary.md)
