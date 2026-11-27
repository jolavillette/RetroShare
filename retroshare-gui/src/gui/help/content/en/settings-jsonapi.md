# JSON API

The JSON API is a local HTTP interface to your node. Another program sends it
requests, each one signed with a token, and can then do what the RetroShare
window does: read and send messages, manage friends, start downloads. The
[web interface](webui.md) uses it, and so can scripts and other clients.

This page of the preferences starts and stops the JSON API server and manages
the tokens. Unless you know what you are doing, you do not need to change
anything here: the web interface registers its own token when you enable it in
the [Webinterface](settings-webui.md) page.

The page exists only when RetroShare was built with the JSON API, which is the
default.

## Layout

The group **JSON API Server** holds, from top to bottom:

- the checkbox **Enable RetroShare JSON API Server** and, on the same line,
  **Status:** with a light that is on while the server runs;
- **Listen Address:** and **Port:**;
- **Token:** with the buttons **Add** and **Remove**;
- the list **Authenticated Tokens:**;
- the list **Registered services:**.

The button **Apply settings** is under the group.

## Buttons

- **Enable RetroShare JSON API Server**: when checked, the server is started
  with the address and port of the page, and the choice is saved. When
  unchecked, the server is stopped, and the address, the port, the token
  buttons, the token list and **Apply settings** are greyed. Stopping the
  server also stops the web interface, which runs on top of it. Off by
  default. The box shows whether the server is running now. When RetroShare
  starts, the server is started only if the web interface is enabled in the
  [Webinterface](settings-webui.md) page; otherwise check the box again after
  each start.
- **Listen Address:** the IPv4 address on which the server accepts
  connections. `127.0.0.1`, the default, means that only programs running on
  this computer can reach it. The new address is used after you click
  **Apply settings**.
- **Port:** the TCP port of the server, from 1024 to 65535; 9092 by default.
  The new port is used after you click **Apply settings**. The address and
  the port are remembered by your node.
- **Token:** field: a token written `user:password`, where the user and the
  password are made of letters and digits only. The field refuses any other
  character. Click a line of **Authenticated Tokens:** to copy it into this
  field.
- **Add**: authorises the token of the field. A client that presents this user
  name and password is then allowed to use the API. Adding a token whose user
  already exists replaces its password.
- **Remove**: revokes the token whose user name is in the field. To remove a
  token of the list, click it, then click **Remove**.
- **Authenticated Tokens:** the tokens currently authorised, with their
  passwords in clear text. The token named `webui` is the one of the web
  interface; its password is the one set in the
  [Webinterface](settings-webui.md) page.
- **Registered services:** the extensions currently served through the API
  server, for instance the web interface. The list is read-only.
- **Apply settings**: restarts the server with the address and the port of the
  page. Clients that are connected are cut and have to connect again. If you
  click it several times in a row, the restart is delayed by about ten seconds.

## Security

- Anyone who holds a token controls your node entirely: they can read your
  messages, add friends and reach your files. Give a token only to programs
  you trust and choose a long password.
- The connection is plain HTTP: it is not encrypted, and the token travels with
  each request.
- Keep **Listen Address:** on `127.0.0.1`. Do not expose the port to your
  local network or to the Internet. To use the API from another machine, go
  through an SSH tunnel or a secured proxy.
- Remove the tokens you no longer use.

## How to

### Let a script use your node

1. Check **Enable RetroShare JSON API Server**.
2. Type a token such as `myscript:` followed by a long password in **Token:**
   and click **Add**.
3. Give the script the address `http://127.0.0.1:9092` and the token.

### Change the port

1. Set **Port:** to the new value.
2. Click **Apply settings**.

## See also

- [Webinterface preferences](settings-webui.md)
- [Web interface](webui.md)
- [Command line options](commandline.md)
