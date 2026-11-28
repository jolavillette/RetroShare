# Adding friends

RetroShare connects only to the nodes you have added as friends. To become
friends, two people exchange their **RetroShare IDs**: each one gives their ID
to the other, and each one adds the ID received. The connection starts when
both sides have done it.

Your own ID is on the [Home](home.md) page. The ID of a friend is added with
the friend wizard, opened by the **Add friend** button of the Home page.

## The two formats of an ID

- The **short ID** is the default. It holds the fingerprint of the profile,
  the ID and name of the node, and its addresses. The key of the profile
  itself is not in it: RetroShare obtains it when the two nodes connect.
- The **old certificate format** is the long text used up to version 0.6.5. It
  also holds the complete public key of the profile, with its signatures if
  asked. Use it with a friend whose RetroShare does not understand the short
  ID. You get it with **Use old certificate format** in the share menu of the
  [Home](home.md) page.

The friend wizard accepts both. Some options of its last page exist only for
the old format, because they need the key of the profile.

## Dialogs

### Connect Friend Wizard, page RetroShare ID

- The text box receives the ID of your friend. Paste the short ID or the
  whole certificate, nothing else. A `retroshare://certificate` link is
  accepted too. Do not paste a PGP key here.
- **Paste**: pastes the content of the clipboard in the box.
- **Open**: loads the ID from a `.rsc` file.
- The line under the box says whether the text is a **Valid Retroshare ID**
  (short format) or a **Valid certificate** (old format), or **Not a valid
  Retroshare certificate!**. Hover it to see what the ID contains.
- **Next** is available once the ID is valid. It leads to the page **Make
  Friend**.

If the ID is your own, or cannot be read, the wizard shows the page **Sorry,
some error appeared** with the reason.

### Connect Friend Wizard, page Make Friend

The page shows what the ID contains under **Peer details**:

- **Name:** the name of the profile.
- **Profile ID:** the ID of the profile key.
- **Location:** the name of the node, followed by the node ID.
- **Signers:** the profiles that have signed this profile key. Only for the
  old certificate format.
- **IP-Addr:** the addresses and ports in the ID, with the dynamic DNS name
  in parentheses. For a hidden node, the onion or I2P address.
- **Known IP:** the address history, when the ID contains one.
- **Key validity:** the validity level of the key. Shown only in advanced
  mode.

The line **This peer is already on your friend list. Adding it might just set
it's ip address.** appears when you already accept this profile. Finishing
the wizard then adds the node to the profile if it is a new one, and updates
its addresses.

**Show Advanced options** shows the options below. The button is absent, and
the options are always visible, when **Enable Advanced Mode** is checked in
the [General preferences](settings-general.md).

- **Add key to keyring**: stores the profile key in your keyring. Unchecking
  it also unchecks **Add as friend to connect with**. It is disabled when the
  key is already in the keyring, and for a short ID, which contains no key.
- **Authenticate friend (Sign PGP Key)**: signs the profile key of your friend
  with your own key. Your profile password is asked. A signature tells your
  other friends that you vouch for this key. It is optional and cannot be
  undone. Not available for a short ID; hidden when you have signed the key
  already.
- **Add as friend to connect with**: accepts connections with this profile.
  It is always on for a short ID and for a profile you already accept.
  Without it, the wizard only stores the key.
- **Add IP to whitelist**: adds the external address of the ID to the list of
  addresses that are always allowed. The list next to it chooses the single
  address or its /24 or /16 range. Not shown when your node or the node of
  your friend is hidden.
- **Add friend to group:** puts the friend in one of your groups. See
  [Network](network.md).
- **Can be used as direct source**: lets RetroShare download directly from
  this friend when they have a file you are downloading. Your friend then
  knows that you download this file. Off by default in the wizard. A warning
  icon tells that the general setting in the
  [Files preferences](settings-files.md) overrides this option.
- **Auto-download recommended files**: downloads automatically the files
  attached to a mail from this friend. Useful between your own nodes.
- **Require whitelist clearance to connect**: refuses a connection from this
  friend when it comes from an address that is not in the whitelist. Refused
  attempts are reported in the [Activity](activity.md) page. Not shown when
  your node or the node of your friend is hidden.

The three last options apply to all the nodes of the profile. You can change
them later, see [Friend details](network-details.md).

**Finish** adds the friend and opens the **Connection Progress** window.

### Connection Progress

The window follows the first connection attempt to the node you have just
added. The lines **Network**, **Connect Status**, **DHT Startup**, **Peer
Lookup** and **UDP Setup** show each step and its result, and the text at the
bottom explains what is happening or what failed.

You can close the window at any time: RetroShare keeps trying to connect in
the background. The same window is opened by **Attempt to connect** in the
[Network](network.md) page.

### Connect Friend Wizard, page RetroShare is better with Friends

This page is opened by **Invite via WebMail** in the share menu of the
[Home](home.md) page. It prepares an invitation message that contains your
own certificate.

- **GMail**, **Yahoo**, **Outlook**, **AOL**, **Yandex**: open the new message
  page of the webmail service in your browser, with the subject **RetroShare
  Invitation** and the invitation text filled in.
- **Email**: opens your e-mail program with the same message.

### Friend request

When a node that is not your friend tries to connect to you, usually because
its owner has added your ID before you added theirs, RetroShare reports it:

- an item appears in the [Activity](activity.md) page, with the buttons
  **Accept Friend Request** and **Peer details**;
- a notification can pop up, depending on the
  [Notify preferences](settings-notify.md).

**Accept Friend Request** opens the friend wizard directly on its last page,
titled **Friend request**. It shows the line **You have a friend request
from** with the name when it is known, the same details and options as the
page **Make Friend**, and the **Accept** button in place of **Finish**.

When RetroShare does not have the profile key of the requester, the request
is handled like a short ID: the options that need the key are disabled.
Accept a request only if you recognise the profile ID or expect this person:
the name alone proves nothing.

An attempt whose certificate has a wrong signature is reported in the
Activity page without the accept button.

### Recommend friends

The entry **Recommend friends to each others** of the share menu of the
[Home](home.md) page opens this dialog. It sends some of your friends the
certificates of other friends, so that they can add each other.

- **Recommend friends**: check the nodes you recommend.
- **To**: check the friends who receive the recommendation.
- **Message:** the text of the mail. A default text is proposed.
- **OK**: sends one mail to each recipient, with the certificate links of the
  recommended nodes. At least one friend must be checked in each list.

The recipient clicks a link in the mail to open the friend wizard on the
recommended node. As always, the two recommended friends must add each other.

To recommend a single node, you can also use **Recommend this node to...** in
the [Network](network.md) page.

## How to

### Add a friend from their ID

1. Open the [Home](home.md) page and click **Add friend**.
2. Paste the ID of your friend and click **Next**.
3. Check the name and the profile ID, then click **Finish**.
4. Send your own ID to your friend if you have not done it yet.

### Add a friend from a certificate link

1. Click a `retroshare://certificate` link in a mail or a chat, or copy it
   and choose **Paste certificate link** in the context menu of the friend
   list of the [Network](network.md) page.
2. The friend wizard opens on the page **Make Friend**. Click **Finish**.

See [RetroShare links](links.md).

## Related settings

- [Network](settings-network.md): your addresses, discovery, hidden service.
- [Notify](settings-notify.md): which connection events are reported.
- [Files](settings-files.md): the general rule for direct downloads.

## See also

- [Getting started](start.md)
- [Home](home.md)
- [Network](network.md)
- [Friend details](network-details.md)
- [Activity](activity.md)
