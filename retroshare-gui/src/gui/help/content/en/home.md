# Home

The Home page is where you find your own RetroShare ID and where you add the
ID of a friend. RetroShare only connects to the nodes you have added as
friends, and a connection needs both sides: you add your friend's ID and your
friend adds yours.

## Layout

- At the top, the line **This is your Retroshare ID. Copy and share with your
  friends!** is followed by your ID. The text is selectable. Hover it to see
  the details of the certificate it contains.
- Below the ID, the block **Did you receive a Retroshare ID from a friend?**
  holds the **Add friend** button.
- The block **Do you need help with Retroshare?** holds the **Open the manual**
  button.

## Buttons

- **?** (help button, top right): shows a short introduction to the page. Click
  the introduction to close it. It is shown automatically the first time you
  open the page.
- **Copy** button, next to the ID (tooltip **Copy your RetroShare ID to
  clipboard**): copies the ID as it is displayed.
- **Share** button (tooltip **Share your RetroShare ID**): opens the share
  menu described below.
- **Add friend**: opens the friend wizard on the page where you paste the ID
  you received. See [Adding friends](start-friends.md).
- **Open the manual**: opens this manual.

## Context menus

### Share menu

The first entries send your ID to someone:

- **Copy your Retroshare ID to Clipboard**: same as the Copy button.
- **Save to File**: writes the ID to a `.rsc` file that you can hand over on
  a USB stick or as an attachment.
- **Send via Email**: opens your e-mail program with a new message that
  contains the ID.
- **Invite via WebMail**: opens the friend wizard on the page that prepares an
  invitation for a webmail service.
- **Recommend friends to each others**: opens a dialog where you choose some of
  your friends and send them the IDs of the others, so that they can connect
  to each other.

The checkable entries change the ID that is displayed, and therefore what the
entries above copy, save or send:

- **Use old certificate format**: shows the long certificate used up to version
  0.6.5 instead of the short ID. Use it for a friend whose RetroShare is too
  old to understand the short format.
- **Include current local IP**: puts the address of your node on your local
  network in the ID. It lets a friend on the same network find you directly.
- **Include current external IP**: puts your public address in the ID.
- **Include my DNS**: puts your dynamic DNS name in the ID, when you have set
  one in the [Network preferences](settings-network.md).
- **Include all IPs history**: puts every address your node has had recently
  in the ID. Off by default.

The four address entries are not shown on a hidden node, which is only
reachable through Tor or I2P and never publishes an IP address.

## How to

### Connect to a friend

1. Click the Copy button and send the ID to your friend by any means: e-mail,
   instant messenger, a file.
2. Ask your friend for their ID.
3. Click **Add friend**, paste the ID and finish the wizard.
4. Wait: the connection starts once both of you have added the other one. Your
   friend then appears in the [Network](network.md) page.

## Related settings

- [Node](settings-node.md): the details of your node and its certificate.
- [Network](settings-network.md): addresses, port, dynamic DNS, hidden
  service.

## See also

- [Getting started](start.md)
- [Adding friends](start-friends.md)
