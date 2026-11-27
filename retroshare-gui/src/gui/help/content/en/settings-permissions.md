# Permissions

RetroShare is made of services: chat, mail, file transfer, forums, channels and
so on. Each service talks separately with each of your friend nodes. This page
of the preferences shows, as a matrix of switches, which service is open with
which friend, and lets you close a service for one friend or for everybody.

A service works between you and a friend only when both sides allow it. When
you switch a service off for a friend, your node stops sending anything of
that service to this friend and ignores what the friend sends for it, and it
tells the friend's node that the service is closed. It works in both
directions: you do not offer the service to that friend and you do not use
theirs.

Every click is applied and saved at once.

## Layout

- At the top, the title **ServicePermissions** and the checkbox
  **hide offline**.
- Below, the matrix:
  - one row per friend node, named "name (location)" and sorted by name. Long
    names are cut after 20 characters;
  - one column per service that runs on your node, with the name of the
    service above it;
  - under each service name, the global switch of the service;
  - at each crossing, the switch of that service for that friend.

Each switch of the matrix has two lights: one for your side and one for your
friend's side. The service passes between you two only when both are on. The
light of a friend who is not connected is shown off, because their node has
not announced its services.

The matrix refreshes every five seconds.

## Buttons

- **hide offline**: when checked, the matrix only shows the friend nodes that
  are connected. Checked each time the page is created.
- A switch at a crossing: click it to deny the service to this friend, click
  again to allow it. The click does nothing while the service is switched off
  globally.
- A global switch: click it to turn the service off for all your friends at
  once; click again to turn it back on for all of them. Turning it back on
  also forgets the friends you had denied one by one.

### Information shown under the mouse

- Over a switch of the matrix, a box shows **Service name:**, **Peer name:**,
  **Peer Id:**, then the friend's side, **Enabled by remote peer** or
  **Disabled by remote peer**, then your side, **Enabled for this peer**,
  **Disabled for this peer** or **Globally switched Off**. The name of the
  friend is underlined in the left column.
- Over a global switch, a box shows the number and the name of the service and
  a short description of what it is for.

## Services

The columns carry the internal names of the services. The usual ones are:

- **disc**: discovery. Friends exchange the addresses and the certificates of
  their common friends. Without it your node relies on the DHT alone to find
  the current address of a friend.
- **chat**: direct chat, distant chat and chat rooms.
- **msg**: mail. Also needed by direct chat between friends.
- **turtle**: anonymous routing. Used by anonymous file transfer and file
  search, distant chat, distant mail and the distant synchronisation of
  channels and the like.
- **heartbeat**: checks that the friend is still alive.
- **ft**: file transfer. Without it you cannot download the files a friend
  shares with you directly; anonymous transfers are not affected.
- **Global Router**: used by distant mail for immediate delivery through
  anonymous tunnels.
- **file_database**: exchange of the lists of shared directories, that is the
  files your friends can browse.
- **serviceinfo**: tells your friends which services are on at your side, and
  the reverse. This service always passes, whatever the switch.
- **bandwidth_ctrl**: speed management.
- **GxsTunnels**: used by distant chat, distant mail and the distant
  synchronisation of channels to carry data through anonymous tunnels.
- **banlist**: exchange of IP filter lists.
- **status**: shares your status (online, away, busy) with friends.
- **gxsid**: exchange of identities. Needed by everything that uses
  identities: chats, forums, mail and so on.
- **gxsforums**: transfer of forums.
- **gxsposted**: transfer of boards.
- **gxschannels**: transfer of channels.
- **gxscircle**: transfer of circles.
- **gxsreputation**: exchange of votes, that is the bans and positive
  opinions on identities.
- **GXS Mails**: used by distant mail for deferred delivery: the mail is
  stored at friends while the recipient is offline.
- **rtt**: measures the round trip time between you and your friends.

A plugin that talks with your friends, such as VOIP, adds its own column.

Be careful: services depend on each other. Turning **turtle** off also stops
all anonymous transfers, distant chat and distant mail. Turning **gxsid** off
breaks everything that needs identities.

## How to

### Stop sharing your file lists with one friend

1. Uncheck **hide offline** if the friend is not connected.
2. Find the row of the friend and the column **file_database**.
3. Click the switch at the crossing. Hover it: the box now reads
   **Disabled for this peer**.

### Turn a service off for everybody

1. Click the global switch under the name of the service.

## See also

- [Network](network.md)
- [Friend details](network-details.md)
- [Network preferences](settings-network.md)
- [Glossary](glossary.md)
