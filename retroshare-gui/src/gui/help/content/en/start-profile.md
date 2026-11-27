# Profiles and nodes

A **profile** is your key pair (a PGP key), protected by a password. It is
your stable identity towards your friends: when a friend adds you, it is your
profile that they accept.

A **node** is one installation of RetroShare that runs with a profile. Each
node has its own name, its own certificate and its own data (friends,
settings, downloaded content). One profile can have several nodes, for
instance one on a desktop and one on a laptop: your friends see them grouped
under your profile name in the [Network](network.md) page.

When RetroShare starts, it shows the login dialog if at least one node exists
on the computer, and the profile creation dialog otherwise.

## Dialogs

### RetroShare (login dialog)

- The list at the top holds the nodes that exist on this computer, each one
  written as profile name, the last eight characters of the profile key ID,
  and node name. The node used last is preselected.
- **Password**: the password of the profile. After a wrong password the field
  is cleared and shows **Wrong password**.
- **Log In**: starts the selected node. The Enter key does the same.
- **New Profile/Node**: opens the profile creation dialog, to create a new
  profile or to add a node to an existing profile. The existing profiles and
  nodes are not changed.

Closing the dialog quits RetroShare.

### Create new profile and new Retroshare node

The dialog creates a profile and its first node in one step. With **Use
existing profile...** checked, its title becomes **Create new Retroshare
node** and it creates only a node.

The basic fields:

- **Node type:** chooses how the node reaches the network. The line above the
  list explains the selected type.
    - **Standard node**: your node connects directly to your friends, who
      therefore see your IP address. It can also connect to hidden nodes when
      Tor runs on your machine.
    - **Hidden node (over Tor)**: all the traffic goes through the Tor
      network and your IP address stays hidden from your friends. RetroShare
      starts and configures Tor by itself. Tor must be installed on the
      computer: without it the creation stops with the message **Tor is not
      available**.
    - **Hidden node (Tor/I2P - Manually configured)**: a hidden node for which
      you set up the Tor hidden service or the I2P tunnel yourself. For
      advanced users.
- **advanced options**: shows the additional fields described below.
- **Username**: the name of the profile, at least 3 characters. Your friends
  see it. Unless you fill **Chat name** in the advanced options, a first
  identity with the same name, signed by the profile, is created too.
- **Password** and **Password again**: the password of the profile, at least
  3 characters, typed twice. It is asked at every start and protects all your
  data. It cannot be recovered if you forget it.
- The progress bar collects randomness for the keys. Move the mouse around
  until it reaches at least 20%.
- **Go!**: creates the profile and the node, then starts RetroShare. The
  button appears only when all the fields are valid and enough randomness has
  been collected; until then a message says what is missing. Creating the
  key can take a while on a slow computer.

A green mark next to a field shows that its content is accepted.

The fields shown by **advanced options**:

- **Use existing profile...**: creates a new node for a profile that already
  exists instead of creating a new profile. The **Username** and **Password
  again** fields disappear, and **Password** must be the password of the
  existing profile.
- The profile list, shown with **Use existing profile...** when the computer
  holds at least one profile: the profile that the new node belongs to.
- **Import profile**, shown with **Use existing profile...**: loads a profile
  from a `.asc` file exported on another computer, and adds it to the profile
  list.
- **Node name**: the name of this installation, at least 3 characters. The
  default is `My computer`. Give different names to the nodes of one profile.
- **Chat name**: the name of the first identity, signed by the profile. Leave
  it empty to create no identity now.
- The key length list: **Default (2048 bits, recommended)**, **High (3072
  bits)** or **Very high (4096 bits)**. It applies to a new profile only.

The fields shown for **Hidden node (Tor/I2P - Manually configured)**:

- **Tor/I2P address**: the onion address of the Tor hidden service, or the
  I2P address, that you have configured for this node. You can leave it empty
  and set it later in the [Network preferences](settings-network.md); the
  node does not work until it is set.
- **Port**: the port of the hidden service, between 1024 and 65535. It can be
  changed later.
- **Use I2P**: lets RetroShare create the I2P tunnel itself through the SAMv3
  interface of your I2P router. The address field then takes the address of
  the I2P router (your own machine when left empty) and **Port** is disabled.

### Setting up Tor...

A node of the type **Hidden node (over Tor)** shows this window at every
start, between the login and the main window, while RetroShare launches Tor
and publishes the hidden service of the node.

Two lines are labelled **Tor status:**.

- The one that shows **Unknown**, **Tor offline** or **Tor ready** is the
  state of Tor itself.
- The one that shows a percentage is the progress of Tor while it joins the
  Tor network. It shows **[Waiting for Tor...]** until Tor answers.
- **Hidden address:** the port of the hidden service, followed by the local
  address and port it forwards to. It shows **[Not ready]** until the service
  exists.
- **Onion address:** the onion address of your node, once it is known.

The window closes by itself when Tor is ready. If Tor cannot be started, an
error message is shown and RetroShare quits.

## How to

### Create a second node for your profile on another computer

1. On the first computer, export the profile to a file from the
   [Node preferences](settings-node.md).
2. Copy the file to the second computer. It is encrypted with your password.
3. On the second computer, start RetroShare, check **advanced options**, then
   **Use existing profile...**, and click **Import profile**.
4. Type the password of the profile, choose a **Node name** and click **Go!**.
5. Each friend must add the new node too: send them the RetroShare ID shown
   in its [Home](home.md) page. Add your own first node the same way to
   connect your two nodes together.

### Add a second profile on this computer

1. In the login dialog, click **New Profile/Node**.
2. Fill the dialog without checking **Use existing profile...**.

## Related settings

- [Node](settings-node.md): the details of your node, the export of the
  profile.
- [Network](settings-network.md): addresses, ports, hidden service
  configuration.

## See also

- [Getting started](start.md)
- [Adding friends](start-friends.md)
- [Command line options](commandline.md)
- [Glossary](glossary.md)
