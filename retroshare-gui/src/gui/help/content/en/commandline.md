# Command line options

RetroShare is normally started from its icon and needs no option. The program
is named `retroshare-gui`. When you start it from a terminal or from a
shortcut, you can add the options below to change where it keeps its data, how
it listens on the network, how it looks, or to hand a link to a RetroShare
that is already running.

## Syntax

`retroshare-gui [options] [links and collection files]`

Each option has a short form (one dash and one letter) and most have a long
form (two dashes and a word). An option that takes a value is followed by the
value: `-p 7812`.

When the value of an option is missing or cannot be read, RetroShare prints
the list of options and stops.

## Options

### Help

- `-h`, `--help`: shows the list of options in a window and in the terminal,
  then stops.

### Data location

- `-c <directory>`, `--base-dir <directory>`: uses this directory instead of
  the default data directory (see below). Profiles, nodes, keys and settings
  are read from it and written to it. Use it to keep a fully separate
  RetroShare, for example on a removable disk or for a test.
- `-R`, `--reset`: resets the settings of the interface (window layout,
  appearance, notification choices...) of the node that is proposed at login
  to their default values. Your profile, friends, messages and files are not
  touched.

### Network

- `-p <port>`, `--port <port>`: listens on this port for this session.
- `-i <address>`, `--ip-address <address>`: uses this local IP address
  (`nnn.nnn.nnn.nnn`) for this session. Useful on a computer with several
  network interfaces.
- `-u`, `--udp`: does not open the listening port for incoming TCP
  connections; the node is reached through UDP only.
- `-o <mode>`, `--opmode <mode>`: starts in the given operating mode: `Full`,
  `NoTurtle`, `Gaming` or `Minimal`. They are the modes **Normal Mode**, **No
  Anon D/L**, **Gaming Mode** and **Low Traffic** of the status bar, see
  [Main window, tray icon and status bar](mainwindow.md). Given while
  RetroShare is already running, the option changes the mode of the running
  RetroShare.
- `-t <path>`, `--tor <path>`: for a hidden node that uses the bundled Tor,
  the full path of the Tor program to run instead of the one RetroShare finds
  by itself.

### Appearance

- `-L <code>`, `--lang <code>`: uses this language for the interface, for
  example `fr` or `de`. With an unknown code RetroShare prints the list of the
  available codes and stops.
- `-g <style>`: uses this interface style, for example `fusion` or `windows`.
  The styles available depend on your system; they are the ones listed in the
  [Appearance preferences](settings-appearance.md).
- `-k <file>`: uses this style sheet.

Without these options RetroShare uses the choices saved in the preferences.

### Logs

- `-s`, `--stderr`: writes the messages of the RetroShare engine to the
  terminal instead of the log file of the node.
- `-d <level>`, `--debug-level <level>`: sets how much the interface logs:
  `debug`, `info`, `notice`, `warn`, `error` or `off`.
- `-l <file>`, `--log-file <file>`: writes the log of the interface to this
  file instead of the terminal.

### Links and files

- `-r <link>`, `--link <link>`: a `retroshare://` link to open.
- `-f <file>`, `--rsfile <file>`: a collection file (`.rscollection`) to
  open.
- Links and collection files can also be given without any option, at the end
  of the command. This is the most reliable form.

Links and files are not opened by the RetroShare you are starting: they are
passed to the RetroShare that is already running, see below.

## Passing a link to a running RetroShare

When RetroShare is already running and you start `retroshare-gui` again with
a link, a collection file or `-o`, the second program does not start. It hands
its arguments to the running RetroShare and stops. The running RetroShare then
opens each link as if you had clicked it, and opens each collection file in
the collection dialog. See [RetroShare links](links.md) for what each kind of
link does.

`retroshare-gui "retroshare://file?name=...&size=...&hash=..."`

`retroshare-gui -o Gaming`

This is what happens when you click a `retroshare://` link in your web
browser once the protocol is registered.

The running RetroShare only listens to such calls when **Accept operating
systems calls** is checked in the
[General preferences](settings-general.md). It is checked by default.

## Where RetroShare keeps its data

Unless you use `-c`, the data directory is:

- Linux and macOS: the directory `.retroshare` in your home directory;
- Windows: the directory `RetroShare` in your application data directory
  (`%APPDATA%\RetroShare`);
- Windows, portable version: the directory `Data` next to the program.
  RetroShare runs as a portable version when a file named `portable` exists
  next to the program.

The data directory contains the PGP keys of your profiles, one directory per
node with its settings, friends and messages, and a directory `extensions6`
for the [plugins](plugins.md) you install yourself. Copy the whole data
directory to back up RetroShare or to move it to another computer. See
[Profiles and nodes](start-profile.md).

## How to

### Run two separate RetroShare installations on one computer

1. Start the second one with its own directory:
   `retroshare-gui -c /path/to/other/directory`.
2. Create or import a profile in it. It shares nothing with the first one.

### Start in a given language once

1. Run `retroshare-gui -L fr`.
2. To keep the language, choose it in the
   [Appearance preferences](settings-appearance.md).

## Related settings

- [General](settings-general.md): starting with the system, starting
  minimised, accepting calls from the operating system, registering the
  `retroshare://` protocol.
- [Network](settings-network.md): the port and addresses used permanently.
- [Appearance](settings-appearance.md): language, style and style sheet.
- [JSON API](settings-jsonapi.md): the port and address of the JSON API.

## See also

- [RetroShare links](links.md)
- [Profiles and nodes](start-profile.md)
- [Web interface](webui.md)
