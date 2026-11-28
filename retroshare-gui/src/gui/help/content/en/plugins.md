# Plugins

A plugin is an optional module, delivered as a separate file, that adds a
feature to RetroShare. RetroShare loads the plugins it finds when it starts.
A plugin can add a page to the main window, a page to the preferences, buttons
to the chat windows, notifications and sounds.

Two plugins come with RetroShare:

- [VOIP](plugins-voip.md): audio and video calls with a friend, from the
  private chat window.
- [FeedReader](plugins-feedreader.md): a reader for RSS and Atom news feeds,
  which can also republish the news in a forum or a board.

Not every RetroShare package contains the plugins. If the
[Plugins preferences](settings-plugins.md) list nothing, your package was
built without them.

## Where RetroShare looks for plugins

At start-up RetroShare searches these directories, in this order:

- on macOS, the `Resources` directory inside the RetroShare application;
- the directory `lib/retroshare/extensions6` located next to the directory
  of the program (for a program in `/usr/bin`, this is
  `/usr/lib/retroshare/extensions6`), then the directory
  `lib/retroshare/extensions6` inside the directory of the program;
- the directory `extensions6` of your RetroShare data directory:
  `~/.retroshare/extensions6` on Linux and macOS,
  `%APPDATA%\RetroShare\extensions6` on Windows, `Data\extensions6` for the
  portable Windows version. RetroShare creates this directory if it does not
  exist. This is where you put a plugin that you add yourself.

A plugin file ends with `.so` on Linux, `.dll` on Windows and `.dylib` on
macOS. When the same plugin is found in two directories, only the first one
is loaded.

The directories that were searched are listed under **Plugin look-up
directories** in the [Plugins preferences](settings-plugins.md).

## Authorising a plugin

A plugin runs inside RetroShare with all its rights, so RetroShare only loads
the plugins you have authorised. It recognises a plugin by the hash
(fingerprint) of its file, and remembers the hashes you accepted.

- The very first time RetroShare runs, the plugins it finds are not loaded.
  Open the [Plugins preferences](settings-plugins.md), click **Enable** on the
  plugins you want and restart RetroShare.
- When a plugin file has changed, or when RetroShare itself has been updated,
  the remembered authorisations no longer apply. At the next start, a dialog
  titled **Unregistered plugin/executable** shows the file and its hash for
  each plugin. Answer **Yes** to authorise the plugin, **No** to refuse it.
  You can change your mind later in the preferences.
- **Disable** in the preferences refuses a plugin. Enabling and disabling take
  effect after a restart.

A plugin must also have been built for your version of RetroShare. One that
was not appears in the preferences with the mark **[loading problem]** and a
status such as **Missing symbol. Wrong version?** or **Loading error.**: get
the plugin that matches your RetroShare.

The option **Authorize all plugins** of the preferences loads every plugin
without checking its hash. It is meant for people who develop plugins. Leave
it off: the hash check protects you from a plugin file that someone has
replaced with a malicious one.

## What a plugin adds

- A page of the main window, with its button among the page buttons:
  FeedReader does.
- A page in the preferences: both **VOIP** and **FeedReader** do. The
  **Configure** button of a plugin in the
  [Plugins preferences](settings-plugins.md) opens it too.
- Buttons in the chat windows: VOIP adds the call buttons to the private chat
  window.
- Toasters, entries in the Activity page and sounds, which you turn on or off
  with the others in the [Notify](settings-notify.md) and
  [Sound](settings-sound.md) preferences.
- A service that your friends' nodes must also run for the feature to work
  between you: for a call, your friend needs the VOIP plugin too.

## How to

### Add a plugin

1. Quit RetroShare.
2. Copy the plugin file into the `extensions6` directory of your data
   directory.
3. Start RetroShare. If it asks about an unregistered plugin, answer **Yes**.
4. Otherwise open [Preferences > Plugins](settings-plugins.md), click
   **Enable** on the plugin and restart RetroShare.

### Remove an old copy of a plugin

After an update of RetroShare, an old plugin file left in your `extensions6`
directory may no longer match the program. Delete it there; the plugins
delivered with the program are found in the program's own directory.

## Related settings

- [Plugins](settings-plugins.md): the list of plugins with their status,
  **Enable**, **Disable**, **Configure**, **About**, and **Authorize all
  plugins**.

## See also

- [VOIP](plugins-voip.md)
- [FeedReader](plugins-feedreader.md)
- [Command line options](commandline.md)
