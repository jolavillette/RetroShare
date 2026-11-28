# Webinterface

The web interface lets you use your node from a web browser instead of the
RetroShare window. See [Web interface](webui.md) for what it offers. This page
of the preferences switches it on, sets its password and opens it in your
browser.

The web interface runs on top of the [JSON API](settings-jsonapi.md): it is
served on the address and port of the JSON API server, and its password is
stored there as the token of the user `webui`. Enabling the web interface
starts the JSON API server if needed.

The page exists only when RetroShare was built with the web interface. This
build option is off by default, so some versions of RetroShare do not have the
page.

## Layout

- At the top, the checkbox **Enable Retroshare WEB Interface** and **Status:**
  with a light that is on while the web interface is served.
- The group **Web parameters** holds **Password:**,
  **Web interface directory:** and the box
  **Allow access from all IP addresses (Default: localhost only)**.
- At the bottom, the buttons **Apply settings** and **Start web browser**, and
  a note.

## Buttons

- **Enable Retroshare WEB Interface**: when checked, the password and the
  directory of the page are applied, the web interface is registered in the
  JSON API server and the server is restarted. When unchecked, the web
  interface is removed from the server, and the group **Web parameters** and
  the two buttons are greyed. The choice is saved and used at the next start.
  Off by default. The box shows unchecked as long as no password has been set.
- **Password:** the password you type in the browser to open the web
  interface. The characters are hidden. The password is taken into account
  when you click **Apply settings** or when you check the box above.
- **Web interface directory:** the directory that holds the files of the web
  interface. The button next to the field (tooltip **Select directory for
  webinterface files (advanced)**) opens a dialog to choose it. By default it
  is `/usr/share/retroshare/webui/`, or `./webui/` on Windows. Change it only
  if the files are installed somewhere else.
- **Allow access from all IP addresses (Default: localhost only)**: this box is
  greyed and has no effect. The address the web interface listens on is the
  **Listen Address:** of the [JSON API](settings-jsonapi.md) page.
- **Apply settings**: saves the page, sets the password and the directory and
  restarts the server.
- **Start web browser**: opens the web interface in your default browser, at
  the address and port of the JSON API server, for instance
  `http://127.0.0.1:9092/index.html`. If the web interface is not enabled, the
  message **The webinterface is not enabled. Enable it in Settings ->
  Webinterface.** is shown instead.

The note at the bottom reminds you that these settings do not concern
`retroshare-service`, the version of RetroShare without a window: it has its
own command line options to activate the web interface and choose the port.

## Security

Do not expose the web interface to the Internet: the connection is not
encrypted and the password protects the whole node. To reach it from another
machine, use an SSH tunnel or a secured proxy. Anyone who knows the password
controls your node.

## How to

### Open the web interface for the first time

1. Type a password in **Password:**.
2. Check **Enable Retroshare WEB Interface**. The light after **Status:**
   turns on.
3. Click **Start web browser**.
4. In the browser, enter the password.

### Change the password

1. Type the new password in **Password:**.
2. Click **Apply settings**.
3. Log in again in the browser with the new password.

## Related settings

- [JSON API](settings-jsonapi.md): the listening address and the port, and the
  list of tokens where the `webui` token appears.

## See also

- [Web interface](webui.md)
- [Command line options](commandline.md)
