# Friend details

Two dialogs describe a friend, because a friend is made of two things:

- the **profile**, the key that identifies the person. The dialog
  **Retroshare profile** shows the key, your trust in it, and the options
  that apply to all the nodes of the friend.
- the **node**, one installation of RetroShare of that person. The dialog
  **Retroshare node details** shows the state of the connection, the
  addresses of the node and its certificate.

Both are opened from the [Network](network.md) page: **Profile details** on a
friend, **Node details** on a node, a double click in the **Keyring** tab.
The **Peer details** button of the items of the [Activity](activity.md) page
opens the node details when the node is known, the profile otherwise.

## Dialogs

### Retroshare node details

The header shows the profile name followed by the node name in parentheses.
**OK** saves the addresses you have changed in the **Connectivity** tab;
**Cancel** closes without saving.

The **Details** tab, under **Node info:**, with the avatar of the node on
the side:

- **Name:** the name and ID of the profile, as a link that opens the
  profile dialog.
- **Node name :** the name of the node.
- **Status :** the state of the connection, with the kind of connection when
  the node is connected.
- **Last Contact :** the last time the node was connected.
- **Retroshare version :** the version of RetroShare that the node runs,
  when it is known.
- **Node ID :** the ID of the certificate of the node.
- **Encryption**: the cipher used by the current connection, or **Not
  connected**.
- **Status message:** the status message your friend has set.

The **Connectivity** tab:

- **Current address:** holds the addresses your node uses to reach this one.
  They are normally filled by the ID of your friend and kept up to date by
  discovery and by the DHT; change them only when you know better.
    - **Local Address** and **Port**: the address of the node on its local
      network. Used when both nodes are on the same network.
    - **External Address** and **Port**: the public address of the node.
    - **Dynamic DNS**: the host name of the node, when your friend uses one.
- For a hidden node, the first line is labelled **Hidden Address** and holds
  the onion or I2P address and its port. The external address and the
  dynamic DNS show **none** and are disabled.
- **List of known addresses:** the addresses at which the node has been seen
  recently.

The third tab shows the ID of the node as text, to copy and hand over to
someone else. Hover the text to see what it contains.

- **short format**: checked by default. The tab is named **Retroshare ID**
  and shows the short ID. Unchecked, the tab is named **Retroshare
  Certificate** and shows the old certificate format, which contains the
  whole profile key. See [Adding friends](start-friends.md).
- **Include signatures**: adds the signatures that other profiles have made
  on the profile key. For the old format only; disabled when the key has no
  such signature.
- **Include IP history**: puts all the known addresses of the node in the ID
  instead of the current ones. Without effect for a hidden node.

### Retroshare profile

**OK** saves the trust level and the content of the **Friend options** tab;
**Cancel** closes without saving. The buttons **Sign this key**, **Accept
connections** and **Deny connections** act at once.

The **Profile info** tab:

- **Name :** the name of the profile.
- **Fingerprint :** the fingerprint of the profile key. To make sure a key
  really belongs to your friend, compare the fingerprint with the one your
  friend reads to you from their own RetroShare.
- **Trust level:** your own opinion of this key: **Unset**, **Unknown**,
  **No trust**, **Marginal**, **Full** or **Ultimate**. It is a personal
  note: it stays on your node and is not used to allow or refuse
  connections.
- **Keysigning:** the **Sign this key** button signs the profile key with
  your own key, after asking your profile password. A signature tells your
  other friends that you vouch for this key. It is optional and cannot be
  undone. The button is disabled once you have signed the key.
- **Accept connections**: shown for a profile that is not your friend. Makes
  your node accept connections with the nodes of this profile. It is better
  to exchange RetroShare IDs, which also carry the addresses needed to
  connect. Disabled when the key type is not supported (only RSA keys are).
- **Deny connections**: shown for a friend. Stops accepting connections with
  this profile: it leaves your friend list at once, without confirmation. The
  key stays in the keyring.
- **This key is signed by :** the profiles that have signed this key. Click
  one to open its own profile dialog.
- **This key has signed your own PGP key**: shown when it is the case.

On your own profile, the trust level and the three buttons are not shown.

The **PGP key** tab shows the public key of the profile as text. This key
identifies all the nodes of the profile. It is not what you exchange to make
friends: for that, use the RetroShare ID of a node.

- **Include signatures**: adds the signatures made by other profiles to the
  text. Disabled when the key has none.

The **Friend options** tab is available only for a profile you accept
connections with. **These options apply to all nodes of the profile:**

- **Use as direct source, when available**: lets RetroShare download
  directly from this friend when they share a file you are downloading. Your
  friend then knows that you download this file. Unchecked, you can still
  download from the file list of the friend on purpose. A warning icon tells
  that the general setting of the [Files preferences](settings-files.md) is
  set to always or never, and overrides this option.
- **Auto-download recommended files from this node**: downloads
  automatically the files attached to a mail whose author is a signed
  identity of this profile. Useful to send files between your own nodes.
- **Require white list clearance**: refuses a connection from this friend
  when it comes from an address that is not in your whitelist. It protects
  against a third party that would relay the connection. A refused attempt
  is reported in the [Activity](activity.md) page, where you can add the
  address to the whitelist.
- **Max upload speed (0=unlimited)**: the highest speed, in kB/s, at which
  your node sends data to this friend.
- **Max download speed (0=unlimited)**: the highest speed, in kB/s, at which
  your node receives data from this friend.

## How to

### Check that a key really belongs to your friend

1. Right-click your friend in the [Network](network.md) page and choose
   **Profile details**.
2. Ask your friend, by phone or face to face, to read the fingerprint of
   their own profile.
3. Compare it with the **Fingerprint :** field.

### Limit the bandwidth used by one friend

1. Open **Profile details** and go to the **Friend options** tab.
2. Set **Max upload speed** and **Max download speed** in kB/s.
3. Click **OK**.

### Fix the address of a friend who does not connect

1. Right-click the node of your friend and choose **Node details**.
2. In the **Connectivity** tab, type the external address and port your
   friend gives you, then click **OK**.
3. Right-click the node again and choose **Attempt to connect**.

## Related settings

- [Permissions](settings-permissions.md): which services each friend may
  use.
- [Files](settings-files.md): the general rule for direct downloads.
- [Network](settings-network.md): your own addresses, the IP filters and the
  whitelist.

## See also

- [Network](network.md)
- [Adding friends](start-friends.md)
- [Activity](activity.md)
- [Glossary](glossary.md)
