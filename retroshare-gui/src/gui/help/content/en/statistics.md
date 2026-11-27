# Statistics

The Statistics window shows what your node is doing on the network: how much
data goes to and comes from each friend, which anonymous tunnels pass through
you, what the routers hold, how fast your friends answer. Nothing in it is
needed for daily use. It helps you understand a slow transfer, a friend that
does not connect or a message that is not delivered yet.

A second, smaller window, the Bandwidth Graph, only draws your total upload and
download speed.

## How to open it

- Right click the RetroShare icon in the system tray and choose
  **Statistics** or **Bandwidth Graph**. On macOS the same entries are in the
  menu of the Dock icon.
- In the status bar of the main window, click the **NAT** or **DHT**
  indicator to open the Statistics window, or the **Down / Up** rates to open
  the Bandwidth Graph. See
  [Main window, tray icon and status bar](mainwindow.md).

The Statistics window is titled "RetroShare Statistics" followed by the name
of your profile and node. Press Esc or close the window to hide it. The curves
keep their history while the window is hidden.

## Layout

The toolbar at the top selects one of the pages described below: **Bandwidth**,
**Turtle Router**, **Identities**, **Global Router**, **Gxs Transport**,
**RTT Statistics** and **DHT**. The **DHT** page does not exist on a hidden
node, which does not use the DHT.

The toolbar follows the button style and icon size chosen in the
[Appearance preferences](settings-appearance.md).

### Graphs

All the curves of the window react to the mouse wheel in the same way:

- wheel alone: stretches or compresses the time axis;
- Shift + wheel: smooths the curve more or less;
- Ctrl + wheel: changes the thickness of the lines.

## Bandwidth

The upper part is a table with one line per connected friend node, and a
**Totals** line. The lower part is a graph of the traffic.

### The table

- **Name**: the name of the friend node.
- **ID**: the identifier of the friend node.
- **In (kB/s)**: the current download speed from this friend.
- **Total In**: the amount of data received from this friend since RetroShare
  started.
- **InMax (kB/s)**: the maximum download speed currently allocated to this
  friend. This value is sent to the friend.
- **InQueue items**: the number of received packets that wait to be processed.
- **Out (kB/s)**: the current upload speed to this friend.
- **Total Out**: the amount of data sent to this friend since RetroShare
  started.
- **OutMax (kB/s)**: the maximum upload speed currently allocated to this
  friend.
- **OutAllowed (kB/s)**: the upload limit requested by the friend, that is the
  friend's own **InMax**.
- **OutQueue items**: the number of packets that wait to be sent to this
  friend.
- **OutQueue bytes**: the total size of the data that waits to be sent.
- **Drain (s)**: the estimated time needed to empty the output queue at the
  current speed. The cell turns orange from 30 seconds and red from 60 seconds:
  the friend receives more slowly than you produce data for them.

The **Totals** line adds up the connected friends. Its **Total In** and
**Total Out** also count the friends that have disconnected since the start.

Click a column header to sort the table. The layout of the columns is
remembered.

### The graph

The graph keeps the last 30 minutes. The controls below it choose what is
drawn:

- **Friend:** **Sum** draws one curve for all the friends together, **All**
  draws one curve per friend, and the name of a friend draws this friend only.
- **Type:** **Up** for the data you send, **Down** for the data you receive.
- **Service:** **Sum** draws all the services together, **All** draws one
  curve per service (chat, file transfer, forums and so on), and the name of a
  service draws this service only. **Friend:** and **Service:** cannot both be
  on **All**: choosing **All** in one sets the other one back to **Sum**.
- **Unit:** **KB/s** draws the speed, **Count** draws the number of packets.
- **Legend:** **Current** shows the latest value next to each curve name. The
  second choice shows the **Average** speed in **KB/s** mode and the **Total**
  number of packets in **Count** mode.
- **Log scale**: uses a logarithmic vertical axis, which makes small values
  visible next to large ones.
- The last list switches the colours of the graph between **Default** and
  **Dark**. The choice also applies to the other graphs of the window the next
  time RetroShare starts.

## Bandwidth Graph

The Bandwidth Graph is a separate small window with the curves of your total
upload and download speed over the last 30 minutes. The text **Since:** gives
the time at which the graph was started or last cleared.

The buttons at the bottom left carry icons only:

- the first button switches the graph between light and dark colours;
- the up arrow and the down arrow show or hide the upload and the download
  curve;
- the eraser button (tooltip **Clears the graph**) empties the graph and
  resets the **Since:** time;
- the settings button (tooltip **Show Settings**, Windows and macOS only)
  shows or hides a panel with a slider (tooltip **Changes the transparency of
  the Bandwidth Graph**, value shown as **% Opaque**) and the **Always on
  Top** check box.

Close the window with Esc on Windows and with Ctrl+W on the other systems.

## Turtle Router

The turtle router is the part of RetroShare that makes anonymous exchanges
possible. A request (a file search, or a tunnel request) is passed from friend
to friend, so that nobody knows where it started. A **tunnel** is the path that
such a request has opened between two nodes that are not friends: the data
then travels through it, relayed by the friends in between, who cannot read it
and do not know the two ends.

### Turtle Router tab

- The graph at the top draws five curves: **TR up** and **TR dn** are the
  tunnel requests you send and receive, **Data up** and **Data dn** are the
  tunnel data you send and receive for yourself, and **Data forward** is the
  data you relay for others.
- **Search requests repartition** and **Tunnel requests repartition** are two
  coloured grids with one line per friend from whom requests arrived. The left
  grid spreads the requests over their age, from 0 to 60 seconds (**Age in
  seconds**). The right grid spreads them over their **Depth**, the number of
  nodes they have already crossed. The number in brackets after the friend's
  name is the **total** for this friend. The colour scale under the grids
  gives the number of requests behind each colour.
- **Turtle router traffic** lists the current speeds: **Tunnel requests Dn**,
  **Tunnel requests Up**, **Incoming file data**, **Outgoing file data** and
  **Forwarded data**.
- **TR Forward probabilities** gives, for each depth (the number in brackets),
  the probability with which your node passes a tunnel request on. RetroShare
  lowers these values when too many requests arrive.

### File transfer tunnels tab

This tab has two sub-tabs.

**Requests** lists what the router has seen recently:

- under **Search requests**, the file searches. **Request ID** is the
  identifier of the search, **From** the friend who passed it to you, **Time**
  its age and **String** the searched words with the number of results found;
- under **Tunnel requests**, the requests to open a tunnel, with the same
  columns. When there are many, only the last 50 are listed.

**Tunnels** lists the open tunnels, grouped by the hash of the file they
serve. The tunnels whose file is not known to your node are grouped under
**Unknown hashes**: they are the tunnels you only relay.

- **Hashes \ Tunnel ID**: the file hash, and under it the identifier of each
  tunnel.
- **Speed**: the current speed in the tunnel.
- **Last Transfer**: the time since data last went through the tunnel.
- **From** and **To**: the two ends of the part of the tunnel that goes
  through your node.

A tunnel identifier is written in italic when nothing has gone through the
tunnel for more than 10 seconds, and in bold when the tunnel is fast.

### GXS sync tunnels tab

Forums, channels and boards can also be synchronised through tunnels with
nodes that are not your friends. The **GXS Groups** list has one line per
group (a forum, a channel...) for which this is in use:

- **Service**: the service the group belongs to.
- **Group ID / Master Key**: the identifier of the group.
- **Policy**: **Passive**, **Active** or **Requesting peers**, when your node
  is looking for nodes that have the group.
- **Status**: **Idle** or **Virtual peers available**, when at least one
  tunnel is open.
- **Last contact**: the time of the last exchange, or **Never**.

Under each group, one line per tunnel shows the identifier of the distant
node, the key of the tunnel, its status (**Tunnel OK** or **Tunnel active**,
followed by **Client** or **Server** according to which side opened it) and
the last contact. The **Random Bias** line under the list is a technical value
of your node.

### Authenticated tunnels tab

Authenticated tunnels are tunnels between two identities, encrypted from end
to end. Distant chat and distant mail use them.

- **Tunnel ID**: the identifier of the tunnel.
- **From**: the identity you use, with its name.
- **To**: the identity at the other end, with its name.
- **Status**: the state of the tunnel, as a number.
- **Total sent** and **Total received**: the amount of data exchanged since
  the tunnel was opened.

## Identities

This page summarises the identities your node knows. The figures are
refreshed every minute.

- **Total identities**: the number of identities known to your node.
- **Usage types**: how many identities were last used for each purpose, for
  example **Message author signature validation** when a forum post signed by
  the identity was checked.
- **Usage per service**: the same count, per service (**Forums**,
  **Channels**, **Chat**...).
- **Identity age (in weeks)**: a histogram of the identities by the date they
  were created, over the last 52 weeks.
- **Last used (hours ago)**: a histogram of the identities by the time they
  were last used.

## Global Router

The global router carries small messages, such as distant mail, to an identity
whose node is not one of your friends. Every node learns which of its friends
is the best next step towards each identity and hands the message over to that
friend.

### Pending packets

The list at the top holds the messages that your node keeps until they are
delivered. Its title **Pending packets** gives their number.

- **ID**: the identifier of the message.
- **Identity Name**: the name of the destination identity.
- **Destinaton**: the identifier of the destination identity.
- **Data status**: where the message stands: Pending, Sent, Receipt OK,
  Ongoing or Done.
- **Tunnel status**: the state of the tunnel used to send it: Unmanaged,
  Pending or Ready.
- **Stored data size**: the size of the message.
- **Data hash**: the fingerprint of the message.
- **Receive time** and **Sending time**: when your node received the message
  and when it last tried to send it.
- **Branching factor**: the number of friends to which the message may be
  handed.
- **Receive time (secs ago)** and **Sending time (secs ago)**: the same two
  times, as a number of seconds.

### Routing matrix

The lower part starts with **Managed keys**, the identities for which your
node accepts messages, each with its service. The **Routing matrix** below has
one line per known identity and one column per friend. The colour of a cell,
read on the scale next to the title, says how good this friend is as a next
step towards this identity. Only a few identities are shown at a time: the
one in the middle is detailed friend by friend, with the value for each
friend in brackets. Turn the mouse wheel over the matrix to move through the
identities.

## Gxs Transport

GXS transport is the way distant mail is delivered when the recipient is not
connected: the encrypted message is published in a shared group, like a post
in a forum that only the recipient can read, and spreads from friend to friend
until it reaches the recipient, who sends a receipt back.

### Pending data items

The list at the top holds the messages you have sent and that are not
acknowledged yet.

- **ID**: the identifier of the message.
- **Destination**: the identity it is sent to.
- **Data status**: the current step, from **Processing** to **Waiting for
  receipt** and **Receipt received**, or a failure (**Receipt signature
  failed**, **Encryption failed**).
- **Data size** and **Data hash**: the size and fingerprint of the message.
- **Sending time**: when the message was sent.
- **Group ID**: the transport group in which it was published.
- **Sending time (secs ago)**: the same time, as a number of seconds.

### Gxs Transport Groups

The lower list holds the transport groups known to your node. Open a group to
see the messages it contains, with their author.

- **Group ID / Author**: the identifier of the group; for a message, its
  author.
- **Publish TS**: the date the message was published.
- **Number of messages**: the number of messages in the group.
- **Local size of data**: the space the group takes on your disk.
- **Subscribed**: **Yes** when your node takes part in the group.
- **Popularity**: the number of your friends that have the group.

### Buttons

- **Refresh** (tooltip **Refresh the transport statistics**): reloads both
  lists. They are also reloaded each time you open the page.

## RTT Statistics

RTT means round trip time: the time a small message takes to go to a friend
and come back. A low and steady value means a good connection.

- **Graph View**: one curve per connected friend, in seconds, over the last
  10 minutes. The legend gives the name of each friend.
- **Table View**: one line per connected friend, with the **Peer Name**, the
  **RTT (ms)** in milliseconds (a question mark as long as no measure is
  available), the **Node ID** and the **IP Address**.

## DHT

The DHT (distributed hash table) is a public directory shared by millions of
computers, the same one BitTorrent uses. RetroShare uses it to find the
current address of your friends and to help two nodes behind routers reach
each other. This page is not shown on a hidden node.

### Your network situation

The lines at the top describe how your node is connected:

- **Net Status:** the overall state, from **NET STATE GOOD!** to warnings such
  as **NET WARNING: Behind NAT** or **NET BAD: Offline**.
- **Network Mode:** whether you are on a local network, behind a router
  (**Behind NAT**) or directly on the Internet (**External IP**).
- **Nat Type:** the kind of router that was detected. A **SYMMETRIC NAT** is
  the hardest to get through.
- **Nat Hole:** how your port is opened in the router: **UPNP FORWARD**,
  **NATPMP FORWARD**, **MANUAL FORWARD** or **NO NAT HOLE**.
- **Connect Mode:** the kinds of connections your node can make.
- **Peer Address:** the address under which the other nodes see you.

The counters next to them summarise your friends: in the DHT, **Unreach:**,
**Online:**, **Offline:** and the total **DHT Peers:**; for the connections,
**Disconnected:**, **Direct:**, **Proxy:** and **Relay:**.

### DHT tab

The list shows the DHT nodes your node knows, friends or not.

- **Bucket**: the group the node is filed in, according to how close its key
  is to yours.
- **IP:Port**: the address of the node.
- **Key**: its identifier in the DHT.
- **Status Flags**: technical flags.
- **Found**: when the node was discovered.
- **Last Sent** and **Last Recv**: when you last wrote to it and heard from
  it.

The **Filter:** field keeps only the lines whose address contains the text
you type.

### Peers tab

One line per friend node, with the state of the search and of the connection
through the DHT.

- **Name** and **PeerId**: the friend and its identifier in the DHT.
- **DHT Status**: **ONLINE**, **Unreachable**, **offline**, **Searching**,
  **Failed** or **Not Active (Maybe Connected!)**.
- **ConnectLogic**: the step of the connection procedure.
- **Connect Status**: **Disconnected**, **Udp Started** or **Connected**.
- **Connect Mode**: **Direct**, **Proxy VIA** or **Relay VIA** followed by the
  node that helps, or **None**.
- **Request Status**: whether a connection request is in progress.
- **Cb Status**: technical state.
- **RsId**: the RetroShare identifier of the friend node.

A **Direct** connection goes straight to your friend. With a **Proxy**, a
third node only helps the two routers to let the connection through. With a
**Relay**, all the data goes through the third node, which is slower.

### Relay tab

The connections that go through a relay: those where you are one end
(**RELAY END**) and those where your node relays for two others (**RELAY
PROXY**). The columns are **Relay Mode**, **Source**, **Proxy**,
**Destination**, **Class**, **Age**, **Last Sent** and **Bandwidth**.

### DHT Graph tab

A curve of the number of RetroShare users found in the DHT over the last
30 minutes.

## Context menus

### DHT list

- **Copy (address) to clipboard**: copies the address of the selected node.
  The entry shows the address itself.

### Pending packets list (Global Router)

- **Details**: opens the details of the destination identity.

### Pending data items list (Gxs Transport)

- **View details**: opens the details of the destination identity.

### Gxs Transport Groups list

- **View details**: on a message, opens the [People](people.md) page on its
  author.

## Related settings

- [Network](settings-network.md): port, DHT, relays, speed limits.
- [Appearance](settings-appearance.md): style of the toolbar; which indicators
  are shown in the status bar.
- [Permissions](settings-permissions.md): which services each friend may use.

## See also

- [Main window, tray icon and status bar](mainwindow.md)
- [Network](network.md)
- [Downloads and uploads](files-transfers.md)
- [Glossary](glossary.md)
