# Plugins

A plugin is a separate library that adds a feature to RetroShare, such as
voice and video calls or a feed reader. This page of the preferences lists the
plugins RetroShare has found, tells whether each one is loaded, and lets you
enable or disable them. See [Plugins](plugins.md) for what the plugins do.

A plugin is loaded when RetroShare starts. Enabling or disabling a plugin
therefore takes effect only after you restart RetroShare.

## Layout

- The group **Plugins** lists the plugin files found, one block per file.
- The group **Plugin look-up directories** lists the directories where
  RetroShare looks for plugin files. The list is read-only.
- The checkbox **Authorize all plugins** is at the bottom.

A block of the list shows the icon of the plugin, its name followed by its
version, and a one-line description. For a file that is not loaded, the name is
the file name followed by **[disabled]** or **[loading problem]**, and the
description is the reason:

- **Plugin disabled. Click the enable button and restart Retroshare**: the
  file has not been accepted yet, or you disabled it.
- **Loading error.**, **Missing symbol. Wrong version?**,
  **No plugin object**, **No API number supplied. Please read plugin
  development manual.**, **No SVN number supplied. Please read plugin
  development manual.**: the file could not be loaded. It is usually a plugin
  built for another version of RetroShare; replace it with the one that comes
  with your version.

## Buttons

Each block has these controls:

- **More** (tooltip **Show more details about this plugin**): shows or hides
  three lines under the block: **File name:** the full path of the plugin
  file, **File hash:** the fingerprint of that file, and **Status:** the
  loading status, for instance **Plugins is loaded.** Hover the status to see
  the error message of a file that failed to load.
- **Enable** (tooltip **Enable this plugin (restart required)**): shown for a
  plugin that is not loaded. It accepts the file, identified by its hash, so
  that it is loaded at the next start. The line "... will be enabled after
  your restart RetroShare." then appears in the block.
- **Disable** (tooltip **Disable this plugin (restart required)**): shown for a
  loaded plugin. It rejects the file: the plugin keeps running until you quit
  and is not loaded at the next start.
- **Configure** (tooltip **Launch configuration panel, if provided by the
  plugin**): opens the settings window of the plugin. Shown only for a loaded
  plugin that has one.
- **About**: opens the information window of the plugin. Shown only for a
  loaded plugin.

Below the list:

- **Authorize all plugins**: shows whether RetroShare loads every plugin file
  without checking its hash. It is meant for people who develop plugins. The
  box is greyed: it cannot be changed from this page.

## How plugins are accepted

For security, a plugin is identified by the hash of its file. A plugin you have
accepted is loaded automatically at every start as long as neither the plugin
file nor the RetroShare program changes. After an update of RetroShare or of
the plugin, the plugin has to be accepted again.

A loaded plugin can add its own pages to the main window, its own page to the
list of the preferences, and its own lines to the
[Notify](settings-notify.md) and [Sound](settings-sound.md) pages.

## How to

### Enable a plugin

1. Find the block of the plugin in the list. Click **More** to check the file
   name if several blocks look alike.
2. Click **Enable**.
3. Quit RetroShare and start it again.

### Disable a plugin

1. Click **Disable** in the block of the plugin.
2. Quit RetroShare and start it again.

## See also

- [Plugins](plugins.md)
- [VOIP](plugins-voip.md)
- [FeedReader](plugins-feedreader.md)
