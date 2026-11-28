# Glossary

The words RetroShare uses for its own notions, in alphabetical order.

## Board

A place where people post links, pictures and short texts, and vote on them.
See [Boards](boards.md).

## Certificate

The long form of a RetroShare ID: a block of text that
contains the public key of a profile and the description of one of its nodes.
The short ID has replaced it in everyday use.

## Channel

A publication with one or a few publishers and many readers. Posts carry
text and files. See [Channels](channels.md).

## Chat room

A discussion among several identities. A room is public or private. Older
texts call it a "lobby". See [Chats](chats.md).

## Circle

A list of identities. A forum, a channel or a board can be restricted to a
circle, so that only its members see it. See [People](people.md).

## Collection

A file with the extension `.rscollection` that lists files and directories
shared in RetroShare. Opening it lets you download the whole set. See
[Files](files.md).

## DHT

Distributed hash table: a public directory, shared with BitTorrent clients,
where your node looks up the current address of your friends' nodes. Only
nodes that already know your ID can find you there. A hidden node does not
use it.

## Distant chat

A private conversation with an identity, not with a friend node. It travels
through tunnels, so neither side learns the address of the other.

## Friend

A node whose RetroShare ID you have accepted. Your node connects only to its
friends. A friendship needs both sides to accept each other. See
[Adding friends](start-friends.md).

## Forum

A threaded discussion. See [Forums](forums.md).

## Hash

The fingerprint of a file, computed from its content. RetroShare identifies a
file by its hash, not by its name, which is why a shared file must be hashed
before others can find it.

## Hidden node

A node that is reachable only through Tor or I2P. It never shows an IP
address to its friends. See [Profiles and nodes](start-profile.md).

## Identity

A name under which you write in chat rooms, forums, channels, boards and
distant mail. An identity is separate from your node: other people see the
identity, not where it comes from. You can own several. A signed identity is
tied to your profile; a pseudonymous one is not. See [People](people.md).

## Node

One installation of RetroShare running under a profile, on one computer. It
used to be called a "location". A profile can have several nodes, for example
one at home and one on a laptop.

## Profile

Your PGP key and the name attached to it. It is what your friends trust when
they add you. It is protected by the password you type at login. See
[Profiles and nodes](start-profile.md).

## Relay

A node that forwards a connection between two friends who cannot reach each
other directly, for example when both are behind a firewall.

## Reputation

The opinion attached to an identity: yours, and the ones of your friends. A
banned identity has its posts hidden and is not forwarded any more. See
[People](people.md).

## RetroShare ID

The text you give to someone so that they can add your node as a friend. See
[Home](home.md).

## RetroShare link

A link that starts with `retroshare://` and points to a file, a forum, a
channel, an identity or another object inside RetroShare. See
[RetroShare links](links.md).

## Tunnel

An anonymous path through a chain of friend connections. Each node on the
path only knows its two neighbours. Tunnels carry file transfers from people
who are not your friends, distant chats and distant mail.

## Turtle router

The part of RetroShare that searches files beyond your friends and that
builds tunnels. See [Statistics](statistics.md).
