# Network

The Network page of the [Preferences](settings.md) sets how your node reaches
the Internet and how other nodes reach it: network mode, port, bandwidth
limits, Tor and I2P proxies, IP filters and relays.

RetroShare only connects to your friends. To connect, two friends must find
each other's address and at least one of them must accept incoming
connections. Most options here deal with one of these two needs.

## Layout

On a standard node the page has four tabs: **Network Configuration**,
**Hidden Service Configuration**, **IP Filters** and **Relay**.

On a hidden node, which is only reachable through Tor or I2P and never shows
an IP address, the **IP Filters** and **Relay** tabs do not exist and several
options of the other tabs are hidden or replaced. Each case is noted below.

## Buttons

### Network Configuration tab

- **Network Mode**: which mechanisms help your friends find your node.
  - **Public: DHT & Discovery**: both mechanisms are on. Recommended on a
    standard node.
  - **Private: Discovery Only**: the DHT is off.
  - **Inverted: DHT Only**: discovery is off.
  - **Dark Net: None**: both are off. Your friends only find you with the
    addresses they already know or that are in your RetroShare ID.

  The DHT is the distributed table of BitTorrent. RetroShare uses it to learn
  the current address of your friends and to answer their connection
  requests; nothing is stored in it. Discovery sends the name and
  identifiers of your friends to the friends you are connected to, to help
  them choose new friends; a friendship is never made automatically.

  On a hidden node the choices are **Discovery On (recommended)** and
  **Discovery Off**; the DHT is always off.
- **NAT**: how incoming connections get through your router. Not shown on a
  hidden node.
  - **Automatic (UPnP)**: RetroShare asks your router to open the port, with
    the UPnP protocol. The port fields are greyed.
  - **Firewalled**: no port is opened. Your friends cannot connect to you
    directly from the Internet; RetroShare relies on outgoing connections
    and on connections set up with the help of the DHT.
  - **Manually Forwarded Port**: you have configured your router yourself to
    forward a port to this computer. Set the port below.
- **Local Address** and **Port:** the address of this computer on your local
  network and the port RetroShare listens on. The address is filled in by
  RetroShare. The port, from 1024 to 65535, can be changed when **NAT** is
  not **Automatic (UPnP)**.
- **External Address** and **Port:** your public address, as RetroShare has
  found it, and the port your friends must use from the Internet. The port
  can be changed when **NAT** is not **Automatic (UPnP)**. On a hidden node
  the field shows **Hidden - See Config** and there is no port.
- **Dynamic DNS**: a host name that always points to your current public
  address, if you have one from a dynamic DNS provider. It can be put in your
  RetroShare ID so that friends find you after your address has changed. Not
  shown on a hidden node.
- **Download limit (KB/s)** and **Upload limit (KB/s)**: the maximum rates
  for the whole program, from 1 to 100000 kB/s, 10000 by default. The
  download total can exceed the limit for a while, for instance when many
  small files are transferred. A very low upload limit may block the low
  priority services such as forums and channels: 50 kB/s is the recommended
  minimum. The line **WARNING: These values don't take into account the
  Relays.** recalls that the relay bandwidth comes on top.

Status indicators, green when the item works (not shown on a hidden node):

- **Local network**: your node has a working local network address.
- **UPnP**: the router has accepted the port opened with UPnP.
- **External ip address finder**: RetroShare has found your public address.

A hidden node shows a single green indicator, **[Hidden mode]**.

Lower part of the tab (not shown on a hidden node):

- **Known / Previous IPs:** the addresses your node has or has had recently,
  marked local or external. They may be given to your friends as contact
  addresses.
- **Clear**: empties this list. Use it when the list holds a wrong or expired
  address that you do not want to pass on to your friends.
- **Allow RetroShare to ask my ip to these DNS servers:** RetroShare asks
  the servers of the list below what your public address is. When the box is
  not checked, RetroShare only learns your address when you connect to
  somebody. Leaving it checked helps when you have few friends, or when you
  are behind a firewall or a VPN. The tooltip of the list shows the external
  addresses that were found.

### Hidden Service Configuration tab, Outgoing Manual Tor/I2P

This sub-tab sets the proxies RetroShare uses to reach hidden nodes. A
standard node can connect to hidden friends as soon as a proxy works.

- **Tor Socks Proxy**: the address and port of the SOCKS proxy of Tor. The
  usual value is 127.0.0.1 and 9050; it must match the configuration of your
  Tor.
- **Tor outgoing Okay**: the indicator is green when something answers on
  that port of your computer. It does not mean that your RetroShare traffic
  goes through Tor: it does only for connections to hidden nodes, or if your
  own node is hidden.
- **I2P Socks Proxy**: the address and port of the SOCKS proxy of I2P. The
  text below explains how to create the client tunnel in the I2P router.
  Both fields are greyed when the automatic I2P of the third sub-tab is
  enabled.
- **I2P outgoing Okay**: same meaning as for Tor.

When Tor is managed by RetroShare itself, the fields are greyed and the tab
says **Tor has been automatically configured by Retroshare. You shouldn't
need to change anything here.**

### Hidden Service Configuration tab, Incoming

Only on a hidden node. This sub-tab holds the address at which your friends
reach you.

- **Local Address** and its port: where the hidden service delivers the
  connections on your computer. The address is 127.0.0.1; the port can be
  changed.
- **Onion Address** (Tor), **I2P Address** (I2P) or **Service Address**: your
  hidden address, ending in `.onion` or `.b32.i2p`, and its port. I2P has no
  port. Tor generates the onion address when you set up the hidden service;
  for I2P, copy the base32 address of the server tunnel you created.
- **Test**: saves the addresses and tries to connect to your hidden address
  through the proxy. The test can produce security warnings about connections
  from 127.0.0.1 in the [Activity](activity.md) feed; they show that the
  connection came through.
- **Tor incoming ok** / **I2P incoming ok**: turns green only after a
  successful test. It means your node can be reached from anywhere through
  the Tor or I2P network.
- **Expected Configuration:** the lines to put in the configuration of Tor,
  or the target of the I2P server tunnel, that match the addresses above.

When Tor is managed by RetroShare itself, the fields are greyed and the test
is not shown.

### Hidden Service Configuration tab, Automatic I2P

Not shown when Tor is managed by RetroShare itself. Here RetroShare drives
I2P by itself through the SAM interface of the I2P router, instead of the
tunnels you create by hand. SAM must be enabled in your I2P router.

- **Enable I2P SAMv3 - changing this requires a restart to fully take
  effect**: switches the feature on. The box is greyed whenever the SAM
  session is not online.
- **advanced mode**: shows the tunnel settings and the server key instead of
  the status text.
- **I2P Instance address**: the address of your I2P router, usually
  127.0.0.1. It is the same value as the address of **I2P Socks Proxy**.
- **SAM accessible**: green when the SAM port of the I2P router (7656)
  answers.
- **Address**: your `.b32.i2p` address. Shown once your node has a server
  key.
- **generate new**: only in advanced mode, once your node has a server key.
  Asks I2P for a new server key, which changes your I2P address.
- The status text (simple mode) says whether RetroShare has set up a client
  or a server tunnel, and in which state SAM is.
- **Tunnel length (in/out)**: the number of I2P routers each incoming and
  outgoing tunnel goes through, from 0 to 7. Longer tunnels are more
  anonymous and slower.
- **Tunnel quantity (in/out)**: the number of tunnels used in parallel, from
  1 to 16.
- **Tunnel variance (in/out)**: a random change of the tunnel length, from
  -1 to 2.
- **Server Key**: the private key of your I2P address, in base64. With a
  key, RetroShare sets up a hidden service for I2P; without one, only a
  client tunnel is created. Paste a key in the box and click **load key** to
  use it. The line below names the algorithms of the key. **load key** and
  **generate new** are greyed whenever the SAM session is not online.
- **Start**, **Restart**, **Stop**: start, restart and stop the SAM session.
  After you change an option, use **Restart** to apply it. The buttons are
  greyed when SAM is not accessible.
- **SAM status**: the indicator shows the state of the session; its tooltip
  gives the detail (offline, establishing a session, online).

### IP Filters tab

Not on a hidden node. IP filtering refuses connections with addresses you do
not trust. It protects you from a node that pretends to be one of your
friends in order to relay and observe the traffic ("masquerading").

- **Activate IP filtering**: switches the filter on. On by default. When it
  is off, the blacklist and its options are greyed.

**IP blacklist** sub-tab: the addresses that are refused. The list is filled
from several sources: masquerading nodes reported by the DHT, ranges you
enter, and ranges reported by your friends. Its columns are **IP Range**,
**Origin** (local, from a friend, from a friend of a friend), **Reason**
and **Comment**. A guess can put the address of a real friend in the list:
in that case move it to the whitelist with the context menu.

- **Ban every IP reported by your friends**: the addresses your friends have
  banned are banned on your node too. This is drastic: it may disconnect you
  from friends and force you to whitelist their addresses. Off by default.
- **Ban every masquerading IP reported by your DHT**: every address the DHT
  reports as masquerading is banned. Also drastic. Off by default.
- **Automatically ban ranges of DHT masquerading IPs starting at** ...
  **IPs**: when this many reported addresses, from 2 to 255, fall in the same
  range, the whole range is banned. On by default, with 2 addresses. Used
  alone, this option protects you well against large scale masquerading.

**IP whitelist** sub-tab: the addresses that are always accepted, even when
they are also in the blacklist. It is filled with the addresses found in the
RetroShare IDs you added by hand, the ranges you enter here and the ones you
accept from a security item of the [Activity](activity.md) feed. For each
friend you can also require that its address be in the whitelist, in the
[Friend details](network-details.md).

**Manual input** group, below both lists:

- The first field takes an IP address, such as `193.190.209.15`. Its
  background changes while the text is not a valid address.
- The number next to it is the size of the range in bits: 32 for the single
  address, 24 or 16 for the whole range that shares its first three or two
  numbers.
- The second field takes a comment of your choice.
- **Add to blacklist**: bans the address or range.
- **Add to whitelist**: always accepts the address or range. This button
  stays available when the filter is off.

Selecting a line of a list copies its address, range and comment to these
fields.

### Relay tab

Not on a hidden node. Two nodes that are both behind a firewall cannot always
connect directly. A relay is a node that passes the traffic between them. It
cannot read that traffic, which is encrypted and authenticated by the two
nodes it serves.

- **Enable Relay Connections**: allows connections that go through a relay,
  and lets your node act as a relay for others within the limits of **Relay
  options**. On by default.
- **Use Relay Servers**: your node also uses the relay servers listed in
  **Relay Server Setup**. It has an effect only when relay connections are
  enabled. Off by default.

**Relay options** group: for each class of nodes, the **Number** of relayed
connections you accept and the **Bandwidth per link** in kB/s. The classes
are **Friends**, **Friends of Friends** and **General** (anyone in the
network). **Total Bandwidth** and the **Total:** line are computed: a link
counts twice its bandwidth, once in each direction. As the page warns, this
bandwidth adds up to the limits of the first tab.

**Relay Server Setup** group:

- The field takes the DHT key of a relay server: 40 hexadecimal characters.
  The box next to it is checked when the key has the right length.
- **Add Server**: adds the key to the list **Server DHT Key**.
- **Remove Server**: removes the selected server.

## Context menus

### A line of the IP blacklist

- **Remove**: removes the line. Only for a range you entered yourself.
- **Move IP ... to whitelist**: removes the line from the blacklist and
  whitelists that single address.
- **Whitelist entire range ...**: same, for the range of 256 addresses
  (`/24`) around it.
- **whitelist entire range ...**: same, for the range of 65536 addresses
  (`/16`) around it.

### A line of the IP whitelist

- **Remove**: removes the line from the whitelist.

## How to

### Accept incoming connections without UPnP

1. In your router, forward a port to this computer.
2. In the **Network Configuration** tab, set **NAT** to **Manually Forwarded
   Port**.
3. Type the port in **Port:** on the **Local Address** line, and the port
   opened on the router in **Port:** on the **External Address** line.

### Connect to a friend whose node is hidden behind Tor

1. Install and start Tor on your computer.
2. Open **Hidden Service Configuration**, sub-tab **Outgoing Manual
   Tor/I2P**.
3. Check the address and port of **Tor Socks Proxy** and wait for the
   indicator **Tor outgoing Okay** to turn green.
4. Add your friend's RetroShare ID as usual. See
   [Adding friends](start-friends.md).

### Unblock a friend that was banned by mistake

1. Open **IP Filters**, sub-tab **IP blacklist**.
2. Right-click the line with the address of your friend.
3. Choose **Move IP ... to whitelist**.

## See also

- [Network](network.md)
- [Friend details](network-details.md)
- [Profiles and nodes](start-profile.md)
- [Node](settings-node.md)
- [Permissions](settings-permissions.md)
- [Glossary](glossary.md)
